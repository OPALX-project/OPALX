# We need to redo the dependency for boost to include python. Because we don't
# want opal main functions to depend on Python.h, we redo the boost find_package 
# here but including boost::python now

# boost::python has the python version hardcoded in the library name (yuck) so
# we resolve the version.

SET(PyBoostLib "python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")
MESSAGE(STATUS "Adding ${PyBoostLib} to boost required list")

set (Boost_LIBS ${Boost_LIBS} ${PyBoostLib})
find_package (Boost 1.66.0
    REQUIRED COMPONENTS ${Boost_LIBS})

# This is a generic "Make a python module" function that can be called to make
# a python module like 
#            python_module(<module_name> <list of source files>)

function(python_module MODULE_NAME PYTHON_SRCS)
    MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
    add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
    set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
    target_link_libraries( ${MODULE_NAME}
        libOPALdynamic
        ${Python3_LIBRARIES}
        ${OPTP_LIBS}
        ${OPTP_LIBRARY}
        ${IPPL_LIBRARY}
        ${GSL_LIBRARY}
        ${GSL_CBLAS_LIBRARY}
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARIES}
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        m
        z
        "-rdynamic"
    )
    target_include_directories(${MODULE_NAME} PUBLIC ${Python3_INCLUDE_DIRS})
    install (TARGETS ${MODULE_NAME} DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
endfunction(python_module)


set(OPAL_SRCS "${OPAL_SRCS}" PARENT_SCOPE)

add_subdirectory(PyCore)
add_subdirectory(PyElements)
