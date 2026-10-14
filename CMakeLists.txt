cmake_minimum_required(VERSION 3.18)
project(gps-sdr-sim C CXX)

option(ENABLE_CUDA "Enable CUDA GPU acceleration" ON)
set(USER_MOTION_SIZE "3000" CACHE STRING "Maximum duration at 10Hz")

set(SRCS gpssim.c l5_code.c)

if(WIN32)
    list(APPEND SRCS getopt.c)
endif()

set(CUDA_SUPPORTED FALSE)
set(USE_NVCC_CUSTOM FALSE)

if(ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CUDA_SUPPORTED TRUE)
    else()
        # Fallback to direct nvcc if available
        find_program(NVCC_EXECUTABLE nvcc HINTS ENV CUDA_PATH PATH_SUFFIXES bin)
        if(NVCC_EXECUTABLE)
            message(STATUS "Found nvcc executable: ${NVCC_EXECUTABLE}")
            set(USE_NVCC_CUSTOM TRUE)
            set(CUDA_SUPPORTED TRUE)
        else()
            message(STATUS "CUDA compiler not found. Building CPU version only.")
        endif()
    endif()
endif()

if(CUDA_SUPPORTED)
    message(STATUS "CUDA acceleration enabled.")
    if(USE_NVCC_CUSTOM)
        set(CUDA_OBJ "${CMAKE_CURRENT_BINARY_DIR}/gpssim_cuda.obj")
        set(NVCC_FLAGS -O3 -arch=native -c "${CMAKE_CURRENT_SOURCE_DIR}/gpssim_cuda.cu" -o "${CUDA_OBJ}")
        if(MSVC)
            list(APPEND NVCC_FLAGS -Xcompiler "/MD /utf-8 /wd4819")
        endif()
        add_custom_command(
            OUTPUT "${CUDA_OBJ}"
            COMMAND "${NVCC_EXECUTABLE}" ${NVCC_FLAGS}
            DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/gpssim_cuda.cu"
                    "${CMAKE_CURRENT_SOURCE_DIR}/gpssim_cuda.cuh"
                    "${CMAKE_CURRENT_SOURCE_DIR}/gpssim.h"
            COMMENT "Compiling gpssim_cuda.cu with nvcc"
        )
        add_executable(gps-sdr-sim ${SRCS} "${CUDA_OBJ}")
        target_compile_definitions(gps-sdr-sim PRIVATE ENABLE_CUDA)
        find_package(CUDAToolkit REQUIRED)
        target_link_libraries(gps-sdr-sim PRIVATE CUDA::cudart)
        if(MSVC)
            target_compile_options(gps-sdr-sim PRIVATE /W3 /O2 /utf-8 /wd4819)
        endif()
    else()
        list(APPEND SRCS gpssim_cuda.cu)
        add_executable(gps-sdr-sim ${SRCS})
        target_compile_definitions(gps-sdr-sim PRIVATE ENABLE_CUDA)
        set_target_properties(gps-sdr-sim PROPERTIES
            CUDA_STANDARD 14
            CUDA_STANDARD_REQUIRED ON
        )
        if(MSVC)
            target_compile_options(gps-sdr-sim PRIVATE
                $<$<COMPILE_LANGUAGE:C>:/W3 /O2 /utf-8 /wd4819>
                $<$<COMPILE_LANGUAGE:CUDA>:-O3>
            )
        endif()
    endif()
else()
    message(STATUS "Building CPU version.")
    add_executable(gps-sdr-sim ${SRCS})
    if(MSVC)
        target_compile_options(gps-sdr-sim PRIVATE /W3 /O2 /utf-8 /wd4819)
    endif()
endif()

if(USER_MOTION_SIZE)
    target_compile_definitions(gps-sdr-sim PRIVATE USER_MOTION_SIZE=${USER_MOTION_SIZE})
endif()

if(UNIX)
    target_link_libraries(gps-sdr-sim PRIVATE m)
endif()

# Unit test for L5 PRN code and Neumann-Hoffman secondary code
add_executable(test_l5_code test_l5_code.c l5_code.c)
if(MSVC)
    target_compile_options(test_l5_code PRIVATE /W3 /O2 /utf-8 /wd4819)
endif()
if(UNIX)
    target_link_libraries(test_l5_code PRIVATE m)
endif()

