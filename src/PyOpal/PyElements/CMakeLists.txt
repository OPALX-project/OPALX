set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}" "elements")
