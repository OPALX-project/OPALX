cmake_minimum_required (VERSION 2.6)

IF(EXISTS ".svn")
    include(FindSubversion)
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} OPAL)
ELSE(EXISTS ".svn")
    set (OPAL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD 
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (OPAL_WC_REVISION \"${GIT_OUT}\")
ENDIF(EXISTS ".svn")

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

MESSAGE("configure: make IpplInfo.h")
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND date
        OUTPUT_VARIABLE DATE_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND uname -a
        OUTPUT_VARIABLE UNAME_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )

set (IPPL_COMPILE_ARCH \"$ENV{IPPL_ARCH}\")
set (IPPL_COMPILE_LINE \"${CMAKE_CXX_FLAGS}\")
set (IPPL_COMPILE_DATE \"${DATE_OUT}\")
set (IPPL_COMPILE_MACHINE \"${UNAME_OUT}\")
set (IPPL_COMPILE_OPTIONS \"${CXXDEFINES}\")
set (IPPL_COMPILE_USER \"adelmann\")
#FIXME: "CVS" version of each file:
#echo `find ${IPPL_ROOT}/src -name '*.cpp' -print` `find ${IPPL_ROOT}/src -name '*.h' -print` | perl ${IPPL_ROOT}/bin/mkver2 >> IpplVersions.h
configure_file(IpplVersions.h.in ${CMAKE_CURRENT_SOURCE_DIR}/IpplVersions.h)

MESSAGE (${CMAKE_CXX_FLAGS})

set (IPPL_SRCS
	DataSource/DataSource.cpp
	DataSource/DataConnect.cpp
	DataSource/DataConnectCreator.cpp
	DataSource/DataSource.cpp
	DataSource/PtclAttribDataSource.cpp
	DataSource/PtclBaseDataSource.cpp
	FieldLayout/FieldLayoutUser.cpp
	FieldLayout/MultiBalancer.cpp
	FieldLayout/VRB.cpp
	Index/Index.cpp
	Meshes/Centering.cpp
	Message/CRC.cpp
	Message/Communicate.cpp
	Message/CommCreator.cpp
	Message/MessageFunctions.cpp
	Message/Formatter.cpp
	Utility/DiscBuffer.cpp
	Utility/DiscConfig.cpp
	Utility/DiscMeta.cpp
	Utility/DiscParticleFunctions.cpp
	Utility/FieldDebugFunctions.cpp
	Utility/ParticleDebugFunctions.cpp
	Utility/Inform.cpp
	Utility/Pool.cpp
	Utility/IpplCounter.cpp
	Utility/IpplStats.cpp
	Utility/IpplTimings.cpp
	Utility/IpplMessageCounter.cpp
	Utility/IpplMemoryUsage.cpp
	Utility/IpplInfo.cpp
	Utility/RNGSimple.cpp
	Utility/RNGXDiv.cpp
	Utility/RNGXCI.cpp
	Utility/RandomNumberGen.cpp
	Utility/Timer.cpp
	Utility/Unique.cpp
	Utility/User.cpp
	Utility/UserList.cpp
	Utility/PAssert.cpp
	Message/CommMPI.cpp
)

set (IPPL_SRCS_FORT
	FFT/fftpack.F
)

set (IPPL_HEADERS_TO_INSTALL
    AppTypes/AntiSymTenzor.h
    AppTypes/AppTypeTraits.h
    AppTypes/dcomplex.h
    AppTypes/GenVektor.h
    AppTypes/SymTenzor.h
    AppTypes/Tenzor.h
    AppTypes/TSVMetaAssign.h
    AppTypes/TSVMetaBinary.h
    AppTypes/TSVMetaCompare.h
    AppTypes/TSVMetaCross.h
    AppTypes/TSVMetaDotDot.h
    AppTypes/TSVMetaDot.h
    AppTypes/TSVMeta.h
    AppTypes/TSVMetaUnary.h
    AppTypes/Vektor.h
	DataSource/DataSource.h
	DataSource/DataConnect.h
	DataSource/DataConnectCreator.h
	DataSource/DataSource.h
	DataSource/PtclAttribDataSource.h
	DataSource/PtclBaseDataSource.h
	FieldLayout/FieldLayoutUser.h
	FieldLayout/MultiBalancer.h
	FieldLayout/VRB.h
    Field/Assign.cpp
    Field/AssignDefs.h
    Field/AssignGeneralBF.cpp
    Field/AssignGeneralIBF.cpp
    Field/Assign.h
    Field/AssignTags.h
    Field/BareField.cpp
    Field/BareField.h
    Field/BareFieldIterator.h
    Field/BCond.cpp
    Field/BCond.h
    Field/BrickExpression.cpp
    Field/BrickExpression.h
    Field/BrickIterator.cpp
    Field/BrickIterator.h
    Field/CompressedBrickIterator.cpp
    Field/CompressedBrickIterator.h
    Field/Field.cpp
    Field/Field.h
    Field/FieldLoc.h
    Field/FieldSpec.h
    Field/GuardCellSizes.cpp
    Field/GuardCellSizes.h
    Field/IndexedBareField.cpp
    Field/IndexedBareField.h
    Field/IndexedField.h
    Field/LField.cpp
    Field/LField.h
    Field/ReductionLoc.cpp
    Field/ReductionLoc.h
	Index/Index.h
    Index/NDIndex.h
	Meshes/Centering.h
	Message/CRC.h
	Message/Communicate.h
	Message/CommCreator.h
	Message/CommMPI.h
    Message/GlobalComm.h
    Message/GlobalComm.cpp
    Message/Message.h
    Message/Message.cpp
    Message/TagMaker.h
    Message/Tags.h
    Message/Formatter.h
    Profile/Profiler.h
    Profile/Profiler.cpp
	Utility/DiscBuffer.h
	Utility/DiscConfig.h
	Utility/DiscMeta.h
	Utility/Inform.h
	Utility/Pool.h
    Utility/Pstring.h
	Utility/IpplCounter.h
	Utility/IpplStats.h
	Utility/IpplTimings.h
	Utility/IpplMemoryUsage.h
	Utility/IpplMessageCounter.h
    Utility/RNGAssignDefs.h
	Utility/RNGSimple.h
    Utility/RNGRand.h
	Utility/RNGXDiv.h
	Utility/RNGXCI.h
	Utility/RandomNumberGen.h
    Utility/SequenceGen.h
	Utility/Timer.h
	Utility/Unique.h
	Utility/User.h
	Utility/UserList.h
	Utility/PAssert.h
	Utility/IpplInfo.h
    PETE/IpplExpressions.h
    PETE/IpplTypeComputations.h
    PETE/PETE.h
    PETE/TypeComputations.h
    Ippl.h
    IpplVersions.h
)

set (MAKE_DEF
    ../test/Makefile.def
)

include_directories (
	${CMAKE_CURRENT_SOURCE_DIR}/DataSource
	${CMAKE_CURRENT_SOURCE_DIR}/FieldLayout
	${CMAKE_CURRENT_SOURCE_DIR}/Index
	${CMAKE_CURRENT_SOURCE_DIR}/Meshes
	${CMAKE_CURRENT_SOURCE_DIR}/Message
	${CMAKE_CURRENT_SOURCE_DIR}/Utility
    ${CMAKE_CURRENT_SOURCE_DIR}/
)

add_library( ippl ${IPPL_SRCS} ${IPPL_SRCS_FORT} )

#foreach(HEADER ${IPPL_HEADERS_TO_INSTALL})
    #string(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
    #install(FILES ${HEADER} DESTINATION include/${DIR})
#endforeach(HEADER)

#FIXME (hack): Instead of copying all headers and cpp files, just create a link
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX}
        OUTPUT_VARIABLE D_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND ln -fs ${CMAKE_SOURCE_DIR}/src ${CMAKE_INSTALL_PREFIX}/include
        OUTPUT_VARIABLE D_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )

install(FILES ${MAKE_DEF} DESTINATION lib)
install(TARGETS ippl DESTINATION lib)
