set (_SRCS
    FFT2D5Algorithm.cpp
    ReferencePath.cpp
    FFT2D5FieldStorage.cpp
    )

set (_PRIVATE_HDRS
    FFT2D5Algorithm.h
    FFT2D5Algorithm.hpp
    ReferencePath.h
    FFT2D5FieldStorage.h
    )

target_sources(opalx PRIVATE ${_SRCS} ${_PRIVATE_HDRS})
