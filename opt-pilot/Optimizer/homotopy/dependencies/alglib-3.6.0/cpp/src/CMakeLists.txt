file( GLOB SOURCES *.cpp )
#message( "Source File are: "${SOURCES} )

add_library (alglib STATIC ${SOURCES} )
