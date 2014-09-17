cmake_minimum_required(VERSION 2.6)

# find the appropriate gtest includes and lib files; find tests from 
# classic_src/*/*.cpp and opal_src/*/*.cpp and fill into the TEST_SRCS variable.

set (GTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include PARENT_SCOPE)
set (GTEST_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gtest/lib PARENT_SCOPE)
set (GTEST_LINK_TARGETS -lgtest -lgtest_main -lpthread PARENT_SCOPE)
set (TEST_EXE opal_unit_tests PARENT_SCOPE)

set(SRC_DIRS classic_src opal_src)

foreach(DIR ${SRC_DIRS})
    file(GLOB TMP_SRC ${DIR}/*/*.cpp)
    message(STATUS "In ${DIR}/*/*.cpp found source files ${TMP_SRC}")
    set(TEST_SRC_LOCAL ${TMP_SRC} ${TEST_SRC_LOCAL})
endforeach(${DIR})

set(TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Main.cpp ${TEST_SRC_LOCAL} PARENT_SCOPE)



