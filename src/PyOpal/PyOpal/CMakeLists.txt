list(APPEND Boost_LIBS python37)
find_package (Boost
    REQUIRED COMPONENTS ${Boost_LIBS})

set (PYUTILS_SRCS Globals.cpp ExceptionTranslation.cpp)

function(python_module MODULE_NAME PYTHON_SRCS)
    MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
    add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
    set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
    target_link_libraries( ${MODULE_NAME}
        OPAL
        ${OPTP_LIBS}
        ${OPTP_LIBRARY}
        ${IPPL_LIBRARY}
        ${GSL_LIBRARY}
        ${GSL_CBLAS_LIBRARY}
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Python3_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        m
        z
        "-rdynamic"
    )
    target_include_directories(${MODULE_NAME} PUBLIC ${Python3_INCLUDE_DIRS})
    install (TARGETS ${MODULE_NAME} DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
endfunction(python_module)

set (POLYNOMIALMAP_SRCS PySquarePolynomialMap.cpp ${PYUTILS_SRCS})
python_module(polynomial_map "${POLYNOMIALMAP_SRCS}")

set (POLYCOEFF_SRCS PyPolynomialCoefficient.cpp ${PYUTILS_SRCS})
python_module(polynomial_coefficient "${POLYCOEFF_SRCS}")

set (PARSER_SRCS PyParser.cpp)
python_module(parser "${PARSER_SRCS}")

set (FIELD_SRCS PyField.cpp ${PYUTILS_SRCS})
python_module(field "${FIELD_SRCS}")

set (BOOSTFIELD_SRCS PyBoostField.cpp ${PYUTILS_SRCS})
python_module(boost_field "${BOOSTFIELD_SRCS}")

set (BINDFIELD_SRCS PyBindField.cpp ${PYUTILS_SRCS})
python_module(bind_field "${BINDFIELD_SRCS}")

set (NDGRID_SRCS PyNDGrid.cpp ${PYUTILS_SRCS})
python_module(ndgrid "${NDGRID_SRCS}")

set (POLYNOMIALPATCH_SRCS PyPolynomialPatch.cpp ${PYUTILS_SRCS})
python_module(polynomial_patch "${POLYNOMIALPATCH_SRCS}")

set (OPALELEMENT_SRCS PyOpalElement.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(opal_element "${OPALELEMENT_SRCS}")

set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}")

set (LINE_SRCS PyLine.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(line "${LINE_SRCS}")

set (TRACKRUN_SRCS PyTrackRun.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(track_run "${TRACKRUN_SRCS}")

set (BEAM_SRCS PyBeam.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(beam "${BEAM_SRCS}")

set (FIELDSOLVER_SRCS PyFieldSolver.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(field_solver "${FIELDSOLVER_SRCS}")

set (TRACKCMD_SRCS PyTrackCmd.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(track "${TRACKCMD_SRCS}")

set (DISTRIBUTION_SRCS PyDistribution.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(distribution "${DISTRIBUTION_SRCS}")

set (RING_DEFINITION_SRCS PyRingDefinition.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(ring_definition "${RING_DEFINITION_SRCS}")

set (LOCAL_CARTESIAN_OFFSET_SRCS PyLocalCartesianOffset.cpp PyOpalObject.cpp ${PYUTILS_SRCS})
python_module(local_cartesian_offset "${LOCAL_CARTESIAN_OFFSET_SRCS}")
