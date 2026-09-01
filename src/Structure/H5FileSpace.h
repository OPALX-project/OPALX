#ifndef H5FILESPACE_H_
#define H5FILESPACE_H_

#include <cstddef>
#include <string>

namespace H5FileSpace {

    /**
     * Rewrite an HDF5 file with only the first @p keepSteps Step# groups.
     *
     * This is a collective operation. Rank zero copies the retained contents to
     * a temporary file and atomically replaces the original after every HDF5
     * handle referring to the original has been closed.
     */
    void compactStepPrefix(const std::string& fileName, std::size_t keepSteps);

}  // namespace H5FileSpace

#endif  // H5FILESPACE_H_
