find_package(Python3 COMPONENTS Interpreter REQUIRED)

execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import h5py, numpy, pandas"
  RESULT_VARIABLE _beambeam_python_dependencies
  ERROR_VARIABLE _beambeam_python_import_error
)
if(NOT _beambeam_python_dependencies EQUAL 0)
  message(
    FATAL_ERROR
      "BeamBeam integration tests require h5py, numpy, and pandas for "
      "${Python3_EXECUTABLE}: ${_beambeam_python_import_error}"
  )
endif()

set(_beambeam_regression_driver
    "${PROJECT_SOURCE_DIR}/sandbox/reduced-order-model/python/run_manufactured_regression.py")
set(_beambeam_regression_work_dir
    "${CMAKE_CURRENT_BINARY_DIR}/manufactured-field-regression")

add_test(
  NAME TestBeamBeamManufacturedFields
  COMMAND
    ${Python3_EXECUTABLE} ${_beambeam_regression_driver}
    --opalx $<TARGET_FILE:opalx_exe>
    --work-dir ${_beambeam_regression_work_dir}
)
set_tests_properties(
  TestBeamBeamManufacturedFields
  PROPERTIES
    LABELS "integration;physics;beambeam"
    PROCESSORS 4
    RUN_SERIAL TRUE
    TIMEOUT 600
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

set(_beambeam_witness_reference_driver
    "${PROJECT_SOURCE_DIR}/sandbox/track12particles/opalx/timed/scan_first_kick_fields.py")
set(_beambeam_witness_reference_work_dir
    "${CMAKE_CURRENT_BINARY_DIR}/witness-reference-offset-regression")

add_test(
  NAME TestBeamBeamWitnessReferenceOffset
  COMMAND
    ${Python3_EXECUTABLE} ${_beambeam_witness_reference_driver}
    --opalx $<TARGET_FILE:opalx_exe>
    --output-dir ${_beambeam_witness_reference_work_dir}
    --case square20_n64_random
    --omp-threads 4
    --force
    --no-plot
    --max-e-magnitude-relative-error 0.08
    --min-e-direction-cosine 0.999
)
set_tests_properties(
  TestBeamBeamWitnessReferenceOffset
  PROPERTIES
    LABELS "integration;physics;beambeam"
    PROCESSORS 4
    RUN_SERIAL TRUE
    TIMEOUT 600
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
