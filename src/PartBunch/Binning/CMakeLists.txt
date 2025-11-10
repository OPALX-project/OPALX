set (_SRCS
    AdaptBins.tpp
    BinHisto.tpp
)

set (HDRS
    AdaptBins.h    
    BinHisto.h
    BinningTools.h
    ParallelReduceTools.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch/Binning")
