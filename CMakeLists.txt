
# set (_HDRS
#     AdaptBins.h
#     AdaptBins.hpp
#     BinningTools.h
#     ParallelReduceTools.h
#     )

# include_DIRECTORIES (
#     ${CMAKE_CURRENT_SOURCE_DIR}
#     )

# add_ippl_headers (${_HDRS})

# install (FILES ${_HDRS} DESTINATION include/binning)

# Add subdirectory for the binning tests
add_subdirectory(test)

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
