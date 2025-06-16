
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
    BinHisto.h
    BinningTools.h
    ParallelReduceTools.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch/Binning")

message (STATUS "Added Binning files.")

