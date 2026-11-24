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

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/PartBunch/Binning")
