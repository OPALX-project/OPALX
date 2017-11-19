IF ( ENABLE_AMR )
#     add_subdirectory(boxlib-amr) still using BoxLib
    add_subdirectory(amrex-only)
    add_subdirectory(ippl-amr)
    add_subdirectory(ippl-bunch)
    IF ( ENABLE_AMR_MG_SOLVER )
        add_subdirectory(trilinos)
    ENDIF ( ENABLE_AMR_MG_SOLVER )
ENDIF (ENABLE_AMR )
