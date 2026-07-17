set (_SRCS
    PicWorkspace.cpp
    )

set (HDRS
    PicWorkspace.h
    PicWorkspace.tpp
    )

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SpaceCharge/Pic")
