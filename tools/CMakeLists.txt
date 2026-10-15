option (ENABLE_BANDRF "Compile BANDRF field conversion scripts" OFF)
if (ENABLE_BANDRF)
    add_subdirectory (BandRF)
endif  ()