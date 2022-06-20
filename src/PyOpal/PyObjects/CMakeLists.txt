# for Parser, we call OpalMain; so must not link to Globals.cpp as we get
# multiple definitions of Opal global objects (i.e. ippl and gmsg)
# hence we don't use ${PYUTILS_SRCS}
set (PARSER_SRCS PyParser.cpp ${PYCORE_SRCS_NO_GLOBALS})
python_module(opal_parser "${PARSER_SRCS}")

set (FIELD_SRCS PyField.cpp ${PYCORE_SRCS})
python_module(field "${FIELD_SRCS}")
