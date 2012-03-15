cmake_minimum_required (VERSION 2.6)

project (DOOM)

set (CMAKE_C_FLAGS
    "-w -O2"
    ${CMAKE_C_FLAGS}
)

set (DOOM_SRCS
    doom.c
)

add_library( doom ${DOOM_SRCS} )

