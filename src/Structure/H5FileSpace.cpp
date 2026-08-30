#include "Structure/H5FileSpace.h"

#include "Utilities/OpalException.h"

#include "Ippl.h"

#include <hdf5.h>
#include <mpi.h>
#include <algorithm>
#include <charconv>
#include <exception>
#include <filesystem>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

namespace {

    struct AttributeCopyContext {
        hid_t destination = H5I_INVALID_HID;
        std::string error;
    };

    void closeAttributeResources(
            hid_t sourceAttribute, hid_t destinationAttribute, hid_t type, hid_t space,
            hid_t creationProperties) {
        if (destinationAttribute >= 0) {
            H5Aclose(destinationAttribute);
        }
        if (creationProperties >= 0) {
            H5Pclose(creationProperties);
        }
        if (space >= 0) {
            H5Sclose(space);
        }
        if (type >= 0) {
            H5Tclose(type);
        }
        if (sourceAttribute >= 0) {
            H5Aclose(sourceAttribute);
        }
    }

    herr_t copyAttribute(hid_t source, const char* name, const H5A_info_t*, void* operationData) {
        auto* context = static_cast<AttributeCopyContext*>(operationData);

        hid_t sourceAttribute      = H5I_INVALID_HID;
        hid_t destinationAttribute = H5I_INVALID_HID;
        hid_t type                 = H5I_INVALID_HID;
        hid_t space                = H5I_INVALID_HID;
        hid_t creationProperties   = H5I_INVALID_HID;

        sourceAttribute = H5Aopen(source, name, H5P_DEFAULT);
        if (sourceAttribute >= 0) {
            type = H5Aget_type(sourceAttribute);
        }
        if (type >= 0) {
            space = H5Aget_space(sourceAttribute);
        }
        if (space >= 0) {
            creationProperties = H5Aget_create_plist(sourceAttribute);
        }
        if (creationProperties >= 0) {
            destinationAttribute = H5Acreate2(
                    context->destination, name, type, space, creationProperties, H5P_DEFAULT);
        }

        const hssize_t points      = space >= 0 ? H5Sget_simple_extent_npoints(space) : -1;
        const std::size_t typeSize = type >= 0 ? H5Tget_size(type) : 0;
        if (destinationAttribute < 0 || points < 0 || typeSize == 0
            || static_cast<unsigned long long>(points)
                       > std::numeric_limits<std::size_t>::max() / typeSize) {
            context->error = "could not prepare root attribute '" + std::string(name) + "'";
            closeAttributeResources(
                    sourceAttribute, destinationAttribute, type, space, creationProperties);
            return -1;
        }

        const std::size_t bytes =
                std::max<std::size_t>(1, static_cast<std::size_t>(points) * typeSize);
        std::vector<unsigned char> value(bytes);
        const bool hasVariableLengthData =
                H5Tdetect_class(type, H5T_VLEN) > 0 || H5Tis_variable_str(type) > 0;

        const herr_t readStatus = H5Aread(sourceAttribute, type, value.data());
        const herr_t writeStatus =
                readStatus >= 0 ? H5Awrite(destinationAttribute, type, value.data()) : -1;
        const herr_t reclaimStatus =
                hasVariableLengthData && readStatus >= 0
                        ? H5Dvlen_reclaim(type, space, H5P_DEFAULT, value.data())
                        : 0;

        closeAttributeResources(
                sourceAttribute, destinationAttribute, type, space, creationProperties);

        if (readStatus < 0 || writeStatus < 0 || reclaimStatus < 0) {
            context->error = "could not copy root attribute '" + std::string(name) + "'";
            return -1;
        }
        return 0;
    }

    bool parseStepIndex(const std::string& name, std::size_t& step) {
        constexpr char prefix[] = "Step#";
        if (name.compare(0, sizeof(prefix) - 1, prefix) != 0 || name.size() == sizeof(prefix) - 1) {
            return false;
        }

        const char* begin = name.data() + sizeof(prefix) - 1;
        const char* end   = name.data() + name.size();
        const auto result = std::from_chars(begin, end, step);
        return result.ec == std::errc() && result.ptr == end;
    }

    std::string linkNameAt(hid_t file, hsize_t index) {
        const ssize_t length = H5Lget_name_by_idx(
                file, "/", H5_INDEX_NAME, H5_ITER_INC, index, nullptr, 0, H5P_DEFAULT);
        if (length < 0) {
            throw OpalException(
                    "H5FileSpace::compactStepPrefix", "could not read an HDF5 root link name");
        }

        std::vector<char> name(static_cast<std::size_t>(length) + 1);
        if (H5Lget_name_by_idx(
                    file, "/", H5_INDEX_NAME, H5_ITER_INC, index, name.data(), name.size(),
                    H5P_DEFAULT)
            < 0) {
            throw OpalException(
                    "H5FileSpace::compactStepPrefix", "could not read an HDF5 root link name");
        }
        return name.data();
    }

    void compactOnRootRank(const std::string& fileName, std::size_t keepSteps) {
        namespace fs = std::filesystem;

        const fs::path original(fileName);
        const fs::path temporary            = original.string() + ".opalx-rewind.tmp";
        const fs::perms originalPermissions = fs::status(original).permissions();
        hid_t source                        = H5I_INVALID_HID;
        hid_t destination                   = H5I_INVALID_HID;

        const auto closeFiles = [&]() {
            if (destination >= 0) {
                H5Fclose(destination);
                destination = H5I_INVALID_HID;
            }
            if (source >= 0) {
                H5Fclose(source);
                source = H5I_INVALID_HID;
            }
        };

        try {
            source = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (source < 0) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not open HDF5 file '" + fileName + "'");
            }

            destination = H5Fcreate(temporary.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            if (destination < 0) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not create temporary HDF5 file '" + temporary.string() + "'");
            }

            AttributeCopyContext attributeContext{destination, {}};
            hsize_t attributeIndex = 0;
            if (H5Aiterate2(
                        source, H5_INDEX_NAME, H5_ITER_INC, &attributeIndex, copyAttribute,
                        &attributeContext)
                < 0) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        attributeContext.error.empty() ? "could not copy HDF5 root attributes"
                                                       : attributeContext.error);
            }

            H5G_info_t rootInfo;
            if (H5Gget_info(source, &rootInfo) < 0) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not inspect HDF5 root objects in '" + fileName + "'");
            }

            for (hsize_t index = 0; index < rootInfo.nlinks; ++index) {
                const std::string name = linkNameAt(source, index);
                std::size_t step       = 0;
                if (parseStepIndex(name, step) && step >= keepSteps) {
                    continue;
                }
                if (H5Ocopy(source, name.c_str(), destination, name.c_str(), H5P_DEFAULT,
                            H5P_DEFAULT)
                    < 0) {
                    throw OpalException(
                            "H5FileSpace::compactStepPrefix",
                            "could not copy HDF5 root object '" + name + "'");
                }
            }

            if (H5Fflush(destination, H5F_SCOPE_GLOBAL) < 0) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not flush temporary HDF5 file '" + temporary.string() + "'");
            }
            if (H5Fclose(destination) < 0) {
                destination = H5I_INVALID_HID;
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not close temporary HDF5 file '" + temporary.string() + "'");
            }
            destination = H5I_INVALID_HID;
            if (H5Fclose(source) < 0) {
                source = H5I_INVALID_HID;
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not close source HDF5 file '" + fileName + "'");
            }
            source = H5I_INVALID_HID;

            std::error_code error;
            fs::permissions(temporary, originalPermissions, error);
            if (error) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not preserve permissions for HDF5 file '" + fileName
                                + "': " + error.message());
            }
            error.clear();
            fs::rename(temporary, original, error);
            if (error) {
                throw OpalException(
                        "H5FileSpace::compactStepPrefix",
                        "could not replace HDF5 file '" + fileName + "': " + error.message());
            }
        } catch (...) {
            closeFiles();
            std::error_code ignored;
            fs::remove(temporary, ignored);
            throw;
        }
    }

}  // namespace

namespace H5FileSpace {

    void compactStepPrefix(const std::string& fileName, std::size_t keepSteps) {
        ippl::Comm->barrier();

        std::string errorMessage;
        if (ippl::Comm->rank() == 0) {
            try {
                compactOnRootRank(fileName, keepSteps);
            } catch (const std::exception& exception) {
                errorMessage = exception.what();
            } catch (...) {
                errorMessage = "unknown error while compacting '" + fileName + "'";
            }
        }

        MPI_Comm comm   = ippl::Comm->getCommunicator();
        int errorLength = static_cast<int>(errorMessage.size());
        if (MPI_Bcast(&errorLength, 1, MPI_INT, 0, comm) != MPI_SUCCESS) {
            throw OpalException(
                    "H5FileSpace::compactStepPrefix", "could not broadcast compaction status");
        }
        if (ippl::Comm->rank() != 0) {
            errorMessage.resize(static_cast<std::size_t>(errorLength));
        }
        if (errorLength > 0
            && MPI_Bcast(errorMessage.data(), errorLength, MPI_CHAR, 0, comm) != MPI_SUCCESS) {
            throw OpalException(
                    "H5FileSpace::compactStepPrefix", "could not broadcast compaction failure");
        }
        if (!errorMessage.empty()) {
            throw OpalException("H5FileSpace::compactStepPrefix", errorMessage);
        }

        ippl::Comm->barrier();
    }

}  // namespace H5FileSpace
