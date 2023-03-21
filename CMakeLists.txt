# ================================================================================
# Project Informationn
#
# ===============================================================================

cmake_minimum_required(VERSION 3.14)
project(apropos VERSION 0.0
                DESCRIPTION "Tapestry: Inference from complex P.f. infections"
                LANGUAGES CXX
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build/cmake/modules/")

# Compilation flags
set(COMPILE_TESTS OFF)


# ================================================================================
# Use C++17
#
# ================================================================================

set(CMAKE_CXX_STANDARD 17)              # required for <filesystem>
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ================================================================================
# Install GoogleTest
#
# ================================================================================

if (COMPILE_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# ================================================================================
# Compile `tapestry`
#
# ================================================================================

# Target executable
set(TARGET tapestry)

# Source files
set(SOURCES 
    src/betabin.cpp
    src/combinatorics.cpp
    src/data.cpp
    src/mcmcs.cpp
    src/models.cpp
    src/parameters.cpp
    src/particle_writers.cpp
    src/particles.cpp
    src/proposals.cpp
    src/sampling.cpp
    src/vcf.cpp
    )

# Add the executable
add_executable(${TARGET} src/main.cpp ${SOURCES})

# --------------------------------------------------------------------------------
# Link htslib
# --------------------------------------------------------------------------------

# Define htslib_FOUND, htslib_INCLUDE_DIRS, htslib_LIBRARIES
find_package(HTSlib 1.4 REQUIRED)

# Link include directories and library to the executable
target_include_directories(${TARGET} PUBLIC ${HTSlib_INCLUDE_DIRS})
target_link_libraries(${TARGET} PUBLIC ${HTSlib_LIBRARIES})


# ================================================================================
# Compile tests
#
# ================================================================================

if (COMPILE_TESTS)
    enable_testing()
    set(TSOURCES 
        tests/test_combinatorics.cpp
        tests/test_particle.cpp
        tests/test_random.cpp
        )
    add_executable(test_tapestry ${TSOURCES} ${SOURCES})

    target_include_directories(test_tapestry PUBLIC ${HTSlib_INCLUDE_DIRS})
    target_link_libraries(test_tapestry GTest::gtest_main ${HTSlib_LIBRARIES})
    
    include(GoogleTest)
    gtest_discover_tests(test_tapestry)
endif()