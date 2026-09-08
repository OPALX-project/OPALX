set (_SRCS
    EndFieldModelManager.cpp
    EndFieldModel.cpp
    Enge.cpp
    Tanh.cpp 
    AsymmetricEnge.cpp
    CompactVector.cpp
    )

set (HDRS
    EndFieldModelManager.h
    EndFieldModel.h
    Enge.h
    Tanh.h
    AsymmetricEnge.h
    CompactVector.h
    )

# Add local include directory
target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Append sources to the main executable target
target_sources(opalx PRIVATE ${_SRCS})

# Optionally install headers
install(FILES ${_HDRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/AbsBeamline/EndFieldModel)
