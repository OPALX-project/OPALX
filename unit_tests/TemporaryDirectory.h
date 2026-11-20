#ifndef OPALX_TEST_TEMPORARY_DIRECTORY_H
#define OPALX_TEST_TEMPORARY_DIRECTORY_H

#include <mpi.h>

#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>

namespace opalx::test {

    /**
     * Own generated test files without sharing filenames with other CTest processes.
     *
     * Atomic directory creation claims a fresh name even for simultaneous invocations
     * of the same executable. Only its creator removes it. MPI fixtures can broadcast
     * one directory to their ranks with create(prefix, comm); they must close shared
     * files and synchronize before reset. Cleanup itself is rank-local and never calls
     * MPI, so destruction during exception unwinding cannot strand peer ranks.
     */
    class TemporaryDirectory {
    public:
        TemporaryDirectory() = default;
        explicit TemporaryDirectory(const std::string& prefix) { create(prefix); }
        ~TemporaryDirectory() { reset(); }

        TemporaryDirectory(const TemporaryDirectory&)            = delete;
        TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

        std::filesystem::path create(const std::string& prefix) {
            return createIn(std::filesystem::temp_directory_path(), prefix);
        }

        std::filesystem::path create(const std::string& prefix, MPI_Comm comm) {
            reset();
            int rank = 0;
            checkMPI(MPI_Comm_rank(comm, &rank));
            std::string error;
            if (rank == 0) {
                try {
                    // CTest's working directory is shared by ranks, unlike node-local /tmp.
                    createIn(std::filesystem::current_path(), prefix);
                } catch (const std::exception& exception) {
                    error = exception.what();
                }
            }
            broadcast(error, rank, comm);
            if (!error.empty()) {
                throw std::runtime_error(error);
            }
            std::string path = directory_m.string();
            broadcast(path, rank, comm);
            directory_m = path;
            return directory_m;
        }

        const std::filesystem::path& path() const { return directory_m; }

        void reset() noexcept {
            if (owner_m) {
                std::error_code ignored;
                std::filesystem::remove_all(directory_m, ignored);
            }
            directory_m.clear();
            owner_m = false;
        }

    private:
        std::filesystem::path createIn(
                const std::filesystem::path& parent, const std::string& prefix) {
            reset();
            const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
            for (unsigned attempt = 0; attempt < 128; ++attempt) {
                const auto suffix    = nextDirectory_m.fetch_add(1, std::memory_order_relaxed);
                const auto candidate = parent
                                       / (prefix + "-" + std::to_string(getpid()) + "-"
                                          + std::to_string(stamp) + "-" + std::to_string(suffix));
                // Existing names are never adopted or removed; mkdir is the ownership claim.
                if (std::filesystem::create_directory(candidate)) {
                    directory_m = candidate;
                    owner_m     = true;
                    return directory_m;
                }
            }
            throw std::runtime_error("could not claim a unique test directory");
        }

        static void checkMPI(int status) {
            if (status != MPI_SUCCESS) {
                throw std::runtime_error("communicating test directory ownership");
            }
        }

        static void broadcast(std::string& value, int rank, MPI_Comm comm) {
            int length = static_cast<int>(value.size());
            checkMPI(MPI_Bcast(&length, 1, MPI_INT, 0, comm));
            if (rank != 0) {
                value.resize(static_cast<std::size_t>(length));
            }
            checkMPI(MPI_Bcast(value.data(), length, MPI_CHAR, 0, comm));
        }

        inline static std::atomic<std::uint64_t> nextDirectory_m{0};
        std::filesystem::path directory_m;
        bool owner_m = false;
    };

}  // namespace opalx::test

#endif
