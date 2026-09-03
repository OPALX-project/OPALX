set (_SRCS
    Pic2d5Solver.cpp
    ReferencePath.cpp
    SliceWorkspace.cpp
    )

set (_PRIVATE_HDRS
    Pic2d5Solver.h
    ReferencePath.h
    SliceWorkspace.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
