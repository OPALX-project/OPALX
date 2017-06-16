include_directories(../dependencies/ )
include_directories(../dependencies/lib )
include_directories(../dependencies/nlopt-2.3/api)
include_directories(../dependencies/alglib-3.6.0/cpp/src)

add_library (Problems STATIC Problems.cpp)
target_link_libraries(Problems wfgProblems alglib)

add_executable (main main.cpp)
target_link_libraries(main nlopt Problems)
