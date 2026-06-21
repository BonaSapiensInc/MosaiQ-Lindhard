cmake_minimum_required(VERSION 3.20)

project(SOPHUS_Lindhard
    VERSION 0.1.0
    DESCRIPTION "C++20 deterministic physics sandbox for finite-temperature Lindhard / RPA"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- strict warnings (build contract from docs/simulator_architecture.md) ---
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# --- release / debug defaults ---
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options($<$<CONFIG:Release>:-O3>)
    add_compile_options($<$<CONFIG:Debug>:-O0>)
    add_compile_options($<$<CONFIG:Debug>:-g>)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_compile_options($<$<CONFIG:Release>:-march=native>)
    endif()
endif()

# --- sophus library (headers + Phase 1 core implementation) ---
add_library(sophus STATIC
    src/core/SincQuadrature.cpp
    src/core/FermiDirac.cpp
    src/engine/Lindhard.cpp
    src/engine/RPA.cpp
    src/engine/PlasmonPoleExtractor.cpp
)
add_library(sophus::sophus ALIAS sophus)

target_include_directories(sophus PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(sophus PUBLIC cxx_std_20)

add_executable(sophus_simulator src/main.cpp)

target_link_libraries(sophus_simulator PRIVATE sophus::sophus)

# --- future unit tests (CTest) ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    enable_testing()
    add_subdirectory(tests)
endif()
