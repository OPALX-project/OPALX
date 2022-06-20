SET(PYUTILS_SRCS ${PYOPAL_ROOT_DIR}/PyCore/Globals.cpp ${PYOPAL_ROOT_DIR}/PyCore/ExceptionTranslation.cpp)

set (_SRCS
    Globals.cpp
    ExceptionTranslation.cpp
    PyOpalObject.cpp
)

set (_HEADERS
    Globals.h
    ExceptionTranslation.h
    PyOpalObject.h
)

set (_SRCS_NO_GLOBALS
    ExceptionTranslation.cpp
    PyOpalObject.cpp
)


include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

# almost every python module needs access to the PyCore sources in order to link things in
SET(PYCORE_SRCS ${OPAL_SOURCE_DIR}/src/PyOpal/PyCore/Globals.cpp 
                ${OPAL_SOURCE_DIR}/src/PyOpal/PyCore/ExceptionTranslation.cpp
                ${OPAL_SOURCE_DIR}/src/PyOpal/PyCore/PyOpalObject.cpp
                parent_scope)
SET(PYCORE_SRCS_NO_GLOBALS
                ${OPAL_SOURCE_DIR}/src/PyOpal/PyCore/ExceptionTranslation.cpp
                ${OPAL_SOURCE_DIR}/src/PyOpal/PyCore/PyOpalObject.cpp
                parent_scope)

add_library(libPyOPALPyCore OBJECT ${_SRCS})
add_library(libPyOPALPyCoreNoGlobals OBJECT ${_SRCS_NO_GLOBALS})

target_include_directories(libPyOPALPyCore PUBLIC ${Python3_INCLUDE_DIRS})
set_property(TARGET libPyOPALPyCore PROPERTY POSITION_INDEPENDENT_CODE 1)

target_include_directories(libPyOPALPyCoreNoGlobals PUBLIC ${Python3_INCLUDE_DIRS})
set_property(TARGET libPyOPALPyCoreNoGlobals PROPERTY POSITION_INDEPENDENT_CODE 1)
