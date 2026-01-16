cmake_minimum_required(VERSION 3.20)

# --- Project Meta Information ---
project(MathLib
    VERSION 1.0.0
    DESCRIPTION "A modern C++20 Math Library"
    LANGUAGES CXX
)

# --- Top Level Check ---
# Check if this project is being built directly or included by another project.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(MATHLIB_IS_TOP_LEVEL ON)
else()
    set(MATHLIB_IS_TOP_LEVEL OFF)
endif()

# --- Options ---
# We use the MATHLIB_IS_TOP_LEVEL variable as the default value.
# You (Top Level) -> Default ON
# User (Subdirectory) -> Default OFF
option(MATHLIB_BUILD_TESTS "Build unit tests" OFF)
option(MATHLIB_BUILD_EXAMPLES "Build example applications" OFF)
option(MATHLIB_BUILD_SANDBOX "Build local sandbox executable" ${MATHLIB_IS_TOP_LEVEL})

# --- C++ Standard ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Library Definition ---
add_library(MathLib INTERFACE)

# --- Include Directories ---
target_include_directories(MathLib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# --- Compiler Warnings (Refined) ---
if(MSVC)
    # MSVC Strict Flags
    target_compile_options(MathLib INTERFACE
        /W4       # Warning level 4
        /permissive- # Standard compliance
    )
else()
    # GCC/Clang Strict Flags (From your list)
    target_compile_options(MathLib INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow            # Warn if a variable shadows another
        -Wconversion        # CRITICAL for Math: Warns on implicit float->int
        -Wsign-conversion   # CRITICAL for Math: Warns on unsigned/signed mix
        -Weffc++            # Effective C++ guidelines
    )
endif()

# --- Subdirectories ---

if(MATHLIB_BUILD_TESTS)
    enable_testing()
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
        add_subdirectory(tests)
    endif()
endif()

if(MATHLIB_BUILD_EXAMPLES)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples")
        add_subdirectory(examples)
    endif()
endif()

# --- Sandbox Executable (Your src/main.cpp) ---
if(MATHLIB_BUILD_SANDBOX AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
    add_executable(Sandbox src/main.cpp)
    # Link the library to the sandbox
    target_link_libraries(Sandbox PRIVATE MathLib)

    # Ensure the sandbox also uses strict warnings (inherited from MathLib)
    # But we set sources specifically for it
    target_include_directories(Sandbox PRIVATE src)
endif()
