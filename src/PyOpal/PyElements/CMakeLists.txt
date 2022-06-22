set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}" "elements")

set (ELEMENT_SRCS PyOpalElement.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(opal_element "${ELEMENT_SRCS}" "elements")

set (RINGDEFINITION_SRCS PyRingDefinition.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(ring_definition "${RINGDEFINITION_SRCS}" "elements")

set (LOCALCARTESIANOFFSET_SRCS PyLocalCartesianOffset.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(local_cartesian_offset "${LOCALCARTESIANOFFSET_SRCS}" "elements")
