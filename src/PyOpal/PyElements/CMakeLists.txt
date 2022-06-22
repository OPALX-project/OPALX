set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}" "elements")

set (ELEMENT_SRCS PyOpalElement.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(opal_element "${ELEMENT_SRCS}" "elements")
