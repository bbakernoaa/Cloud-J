if(WIN32)
  add_test(NAME GenerateOutput
    COMMAND ${CMAKE_COMMAND} -E env cmd /c "$<TARGET_FILE:cloudj_standalone> > actual_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )

  add_test(NAME CompareOutput
    COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol actual_output.txt "${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )
  set_tests_properties(CompareOutput PROPERTIES DEPENDS GenerateOutput)

  if(CLOUDJ_CPP_PORT)
    add_test(NAME GenerateOutputCpp
      COMMAND ${CMAKE_COMMAND} -E env cmd /c "$<TARGET_FILE:cloudj_standalone_cpp> --use-headers-only > cpp_actual_output.txt"
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
    )
    add_test(NAME CompareOutputCpp
      COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol cpp_actual_output.txt "${CMAKE_SOURCE_DIR}/test/expected_output/cpp_reference_output.txt"
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
    )
    set_tests_properties(CompareOutputCpp PROPERTIES DEPENDS GenerateOutputCpp)
  endif()
else()
  add_test(NAME CompareOutput
    COMMAND ${CMAKE_COMMAND} -E env sh -c "$<TARGET_FILE:cloudj_standalone> > actual_output.txt; ${CMAKE_COMMAND} -E compare_files --ignore-eol actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/reference_output.txt"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )

  if(CLOUDJ_CPP_PORT)
    add_test(NAME CompareOutputCpp
      COMMAND ${CMAKE_COMMAND} -E env sh -c "$<TARGET_FILE:cloudj_standalone_cpp> --use-headers-only > cpp_actual_output.txt; ${CMAKE_COMMAND} -E compare_files --ignore-eol cpp_actual_output.txt ${CMAKE_SOURCE_DIR}/test/expected_output/cpp_reference_output.txt"
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
    )
  endif()
endif()

if(CLOUDJ_CPP_PORT)
  add_test(NAME FuzzPropertyParity
      COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/tools/run_property_fuzzer.py --count 10000
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin
  )
endif()

if(CLOUDJ_CPP_PORT)
  add_executable(test_cpp_math
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/cpp/test_cpp_math.cpp
  )
  target_compile_definitions(test_cpp_math PRIVATE CLOUDJ_GPU_MODE)
  set_target_properties(test_cpp_math PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )

  add_executable(test_fortran_math
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/test_fortran_math.F90
  )
  target_link_libraries(test_fortran_math PRIVATE CloudJ_Core)
  set_target_properties(test_fortran_math PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )
endif()



