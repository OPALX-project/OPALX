set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp ${PYUTILS_SRCS} $<TARGET_OBJECTS:libPyOPALPyCoreNoGlobals>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}" "elements")
