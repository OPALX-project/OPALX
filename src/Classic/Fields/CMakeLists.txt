set (_SRCS
  BDipoleField.cpp
  BMultipoleField.cpp
  ConstBField.cpp
  ConstBzField.cpp
  ConstEField.cpp
  ConstEzField.cpp
  EDipoleField.cpp
  EMField.cpp
  NullField.cpp
  StaticElectricField.cpp
  StaticMagneticField.cpp
  Fieldmap.cpp
  FM3DH5Block.cpp
  FM3DH5Block_nonscale.cpp
  FM3DMagnetoStaticH5Block.cpp
  FM3DDynamic.cpp
  FM2DElectroStatic.cpp
  FM2DMagnetoStatic.cpp
  FM2DDynamic.cpp
  FM1DDynamic.cpp
  FM1DDynamic_fast.cpp
  Astra1D_fast.cpp
  Astra1DDynamic.cpp
  Astra1DDynamic_fast.cpp
  FM1DElectroStatic.cpp
  FM1DElectroStatic_fast.cpp
  Astra1DElectroStatic.cpp
  Astra1DElectroStatic_fast.cpp
  FM1DMagnetoStatic.cpp
  FM1DMagnetoStatic_fast.cpp
  Astra1DMagnetoStatic.cpp
  Astra1DMagnetoStatic_fast.cpp
  FM1DProfile1.cpp
  FM1DProfile2.cpp
  FMDummy.cpp
  SectorMagneticFieldMap.cpp
  SectorField.cpp
  )

include_directories (
  Interpolation
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_subdirectory (Interpolation)

add_sources(${_SRCS})