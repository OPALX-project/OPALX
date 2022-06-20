set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp ${PYUTILS_SRCS} $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}")
