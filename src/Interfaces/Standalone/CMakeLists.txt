# Cloud-J/src/CMakeLists.txt

#if("cloudj_standalone" IN_LIST CLOUDJ_EXE_TARGETS)
        add_executable(cloudj_standalone
                cldj_standalone.F90
        )
        target_link_libraries(cloudj_standalone
                PUBLIC CloudJ_Core
        )
        set_target_properties(cloudj_standalone PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
        install(TARGETS cloudj_standalone
                RUNTIME DESTINATION ${RUNDIR}
        )
        if(WIN32)
                target_link_libraries(cloudj_standalone
                        PRIVATE "-Wl,--stack,10485760"
                )
        endif()
#endif()

if(CLOUDJ_CPP_PORT)
        add_executable(cloudj_standalone_cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/standalone/cloudj_standalone.cpp
        )
        set_target_properties(cloudj_standalone_cpp PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
        )
        
        option(CLOUDJ_GPU_MODE_BENCH "Compile C++ standalone in GPU Hermite spline mode" OFF)
        if(CLOUDJ_GPU_MODE_BENCH)
                target_compile_definitions(cloudj_standalone_cpp PRIVATE CLOUDJ_GPU_MODE)
        endif()

        if(WIN32)
                target_link_libraries(cloudj_standalone_cpp
                        PRIVATE "-Wl,--stack,10485760"
                )
        endif()

        if(CLOUDJ_USE_KOKKOS)
                target_link_libraries(cloudj_standalone_cpp PUBLIC Kokkos::kokkos Kokkos::kokkoskernels)
        endif()
endif()

