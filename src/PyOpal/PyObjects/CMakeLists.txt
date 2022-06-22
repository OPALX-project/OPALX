set (PARSER_SRCS PyParser.cpp $<TARGET_OBJECTS:libPyOPALPyCoreNoGlobals>)
python_module(parser "${PARSER_SRCS}" objects)

set (FIELD_SRCS PyField.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(field "${FIELD_SRCS}" objects)

set (LINE_SRCS PyLine.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(line "${LINE_SRCS}" objects)

