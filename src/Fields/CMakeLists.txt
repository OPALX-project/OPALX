set (_SRCS
    # Default Field Classes - DO NOT REMOVE
    BDipoleField.cpp
    BMultipoleField.cpp
    ConstBField.cpp
    ConstBzField.cpp
    ConstEzField.cpp
    EMField.cpp
    NullField.cpp
    StaticElectricField.cpp
    StaticMagneticField.cpp

    # Fieldmaps
    Fieldmap.cpp
    FM2DMagnetoStatic.cpp

    # To be ported from OPAL
    #FM3DH5BlockBase.cpp
    #FM3DH5Block.cpp
    #FM3DH5Block_nonscale.cpp
    #FM3DMagnetoStaticH5Block.cpp
    #FM3DMagnetoStatic.cpp
    #FM3DMagnetoStaticExtended.cpp
    #FM3DDynamic.cpp
    #FM2DElectroStatic.cpp
    #FM2DDynamic.cpp
    #FM1DDynamic.cpp
    #FM1DDynamic_fast.cpp
    #Astra1D_fast.cpp
    #Astra1DDynamic.cpp
    #Astra1DDynamic_fast.cpp
    #FM1DElectroStatic.cpp
    #FM1DElectroStatic_fast.cpp
    #Astra1DElectroStatic.cpp
    #Astra1DElectroStatic_fast.cpp
    #FM1DMagnetoStatic.cpp
    #FM1DMagnetoStatic_fast.cpp
    #Astra1DMagnetoStatic.cpp
    #Astra1DMagnetoStatic_fast.cpp
    #FM1DProfile1.cpp
    #FM1DProfile2.cpp
    #FMDummy.cpp
    #SectorMagneticFieldMap.cpp
    #SectorField.cpp
    )

add_subdirectory (Interpolation)


set (HDRS
    # Default Field Classes - DO NOT REMOVE
    BDipoleField.h
    BMultipoleField.h
    ConstBField.h
    ConstBzField.h
    ConstEzField.h
    EMField.h
    StaticElectricField.h
    StaticMagneticField.h
    NullField.h
    AcceleratingField.h
    OscillatingField.h
    
    # Fieldmaps
    Fieldmap.h
    Fieldmap.hpp
    FM2DMagnetoStatic.h
    
    # To be ported from OPAL
    #FM3DDynamic.h
    #FM3DH5BlockBase.h
    #FM3DH5Block.h
    #FM3DH5Block_nonscale.h
    #FM3DMagnetoStaticH5Block.h
    #FM3DMagnetoStatic.h
    #FM3DMagnetoStaticExtended.h
    #FMDummy.h
    #SectorField.h
    #SectorMagneticFieldMap.h
    #Astra1DDynamic_fast.h
    #Astra1DDynamic.h
    #Astra1DElectroStatic_fast.h
    #Astra1DElectroStatic.h
    #Astra1D_fast.h
    #Astra1DMagnetoStatic_fast.h
    #Astra1DMagnetoStatic.h
    #BSingleMultipoleField.h
    #FM1DDynamic_fast.h
    #FM1DDynamic.h
    #FM1DElectroStatic_fast.h
    #FM1DElectroStatic.h
    #FM1DMagnetoStatic_fast.h
    #FM1DMagnetoStatic.h
    #FM1DProfile1.h
    #FM1DProfile2.h
    #FM2DDynamic.h
    #FM2DElectroStatic.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Fields")