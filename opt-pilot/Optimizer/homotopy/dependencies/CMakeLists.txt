add_subdirectory( nlopt-2.3 )
#add_subdirectory( gtest-1.6.0 )

#if( IS_DIRECTORY lib/boost/ )
	add_subdirectory( lib )
#endif( IS_DIRECTORY lib/boost/ )

add_subdirectory( alglib-3.6.0 )

add_subdirectory( wfgProblems )
add_subdirectory( wfgHypervolume )
