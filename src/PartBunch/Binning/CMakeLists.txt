
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


# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:


#set (_SRCS
#    AdaptBins.tpp
#    BinHisto.tpp
#)

#include_directories (
#    ${CMAKE_CURRENT_SOURCE_DIR}
#    )

#add_opal_sources (${_SRCS})

#set (HDRS
#    AdaptBins.h    
#    # AdaptBins.hpp
#    BinHisto.h
#    # BinHisto.hpp
#    BinningTools.h
#    ParallelReduceTools.h
#)

#install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch/Binning")

#message (STATUS "Added Binning files.")

# add_subdirectory(test)

# Enough in ippl: 
#add_subdirectory(test)
#message (STATUS "Added Binning files.")

set (_SRCS
    AdaptBins.tpp
    BinHisto.tpp
)

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    AdaptBins.h    
    # AdaptBins.hpp
    BinHisto.h
    # BinHisto.hpp
    BinningTools.h
    ParallelReduceTools.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch/Binning")

message (STATUS "Added Binning files.")

