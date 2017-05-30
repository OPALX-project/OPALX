IF ( ENABLE_AMR_SOLVER )
#     add_subdirectory(boxlib-amr) # not AMReX ready, still for BoxLib
    add_subdirectory(boxlib-only)
    add_subdirectory(ippl-amr)
    add_subdirectory(ippl-bunch)
ENDIF (ENABLE_AMR_SOLVER )
