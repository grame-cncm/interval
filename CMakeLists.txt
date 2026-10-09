cmake_minimum_required(VERSION 3.13)

# set the project name
project(TestInterval LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

if (NOTIDY)
else ()
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=-*,clang-*,cppcoreguidelines-*,-cppcoreguidelines-avoid-magic-numbers,performance-*,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,llvm-*,modernize-*,-modernize-use-trailing-return-type,cert-*")
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    if(UNIX)
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
    endif()
endif()

add_definitions(-D_USE_MATH_DEFINES)

file(GLOB INTERVAL_SOURCES CONFIGURE_DEPENDS interval/*.cpp)

add_library(interval STATIC ${INTERVAL_SOURCES})
target_include_directories(interval PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/FaustAlgebra)
# The hosts of the library (the Faust compiler among them) build it as C++17 ; the
# tests may use C++20, the library itself must not (cxx_compat.hh).
set_target_properties(interval PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

# Production code has no multiprecision dependency. MPFR is an explicit opt-in
# test oracle; a normal or cross-compiled build never searches for host packages.
option(INTERVAL_ENABLE_MPFR_TESTS "Check bounds against an independent MPFR oracle" OFF)

# The interval algebra models signed 32-bit integer overflow as two's-complement
# wrapping. GCC and Clang must compile the corresponding signed arithmetic with
# that same contract.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(interval PUBLIC -fwrapv)
    # NaN/inf handling and the steps composing an enclosure must survive user
    # optimization flags. These options also protect header-based affine helpers.
    target_compile_options(interval PUBLIC -fno-fast-math -ffp-contract=off)
    # GCC can otherwise link a fast-math startup object that enables FTZ/DAZ,
    # despite strict compilation of the interval code itself.
    target_link_options(interval PUBLIC -fno-fast-math)
elseif(MSVC)
    target_compile_options(interval PUBLIC /fp:strict)
endif()

# Experimental and numerical analysis executable.
add_executable(TestInterval main.cpp)
target_link_libraries(TestInterval PRIVATE interval)

# Fast deterministic regression tests used by CTest and CI.
add_executable(IntervalTests tests/interval_tests.cpp)
target_link_libraries(IntervalTests PRIVATE interval)

# Boundary regressions remain available without external numerical libraries.
add_executable(DirectedBoundsTests tests/directed_bounds_tests.cpp)
target_link_libraries(DirectedBoundsTests PRIVATE interval)

# The former counterexamples are ordinary inclusion regressions now.
add_executable(FloatConservativenessTests tests/float_conservativeness_tests.cpp)
target_link_libraries(FloatConservativenessTests PRIVATE interval)
add_executable(FloatBoundsTests tests/float_bounds_tests.cpp)
target_link_libraries(FloatBoundsTests PRIVATE interval)

# Integration regressions also compile the affine/public headers as their C++17
# Faust consumers do; the precision loop checks the same typing contract twice.
add_executable(FaustIntegrationTests tests/faust_integration_tests.cpp)
set_target_properties(FaustIntegrationTests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_link_libraries(FaustIntegrationTests PRIVATE interval)

# Validity participates in both domains and their fixpoint ordering. C++17 is
# retained to exercise the same public API used by the Faust compiler.
add_executable(InvalidityTests tests/invalidity_tests.cpp)
set_target_properties(InvalidityTests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_link_libraries(InvalidityTests PRIVATE interval)

# Delay zero must preserve the complete affine form, including empty states;
# the same C++17 regressions cover initialization when a positive delay is possible.
add_executable(DelayTests tests/delay_tests.cpp)
set_target_properties(DelayTests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_link_libraries(DelayTests PRIVATE interval)

enable_testing()
add_test(NAME interval_unit_tests COMMAND IntervalTests)
add_test(NAME directed_bounds_tests COMMAND DirectedBoundsTests)
add_test(NAME float_inclusion_tests COMMAND FloatConservativenessTests)
set_tests_properties(float_inclusion_tests PROPERTIES LABELS "float-inclusion")
add_test(NAME float_bounds_tests COMMAND FloatBoundsTests)
set_tests_properties(float_bounds_tests PROPERTIES LABELS "float-inclusion")
add_test(NAME faust_integration_tests COMMAND FaustIntegrationTests)
add_test(NAME invalidity_tests COMMAND InvalidityTests)
add_test(NAME delay_tests COMMAND DelayTests)

if(INTERVAL_ENABLE_MPFR_TESTS)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(PC_MPFR QUIET mpfr)
    endif()
    find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${PC_MPFR_INCLUDE_DIRS})
    find_path(GMP_INCLUDE_DIR gmp.h HINTS ${PC_MPFR_INCLUDE_DIRS})
    find_library(MPFR_LIBRARY NAMES mpfr libmpfr HINTS ${PC_MPFR_LIBRARY_DIRS})
    find_library(GMP_LIBRARY NAMES gmp libgmp HINTS ${PC_MPFR_LIBRARY_DIRS})
    if(NOT MPFR_INCLUDE_DIR OR NOT GMP_INCLUDE_DIR OR NOT MPFR_LIBRARY OR NOT GMP_LIBRARY)
        message(FATAL_ERROR "INTERVAL_ENABLE_MPFR_TESTS requires MPFR/GMP development files. Set CMAKE_PREFIX_PATH or disable this option.")
    endif()
    add_executable(MPFROracleTests tests/mpfr_oracle_tests.cpp)
    target_include_directories(MPFROracleTests PRIVATE ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    target_link_libraries(MPFROracleTests PRIVATE interval ${MPFR_LIBRARY} ${GMP_LIBRARY})
    add_test(NAME mpfr_oracle_tests COMMAND MPFROracleTests)

    # Optional independent binary32 witnesses, rounding at each target operation.
    add_executable(FloatConservativenessOracleTests tests/float_conservativeness_tests.cpp)
    target_compile_definitions(FloatConservativenessOracleTests PRIVATE INTERVAL_FLOAT_MPFR_ORACLE=1)
    target_include_directories(FloatConservativenessOracleTests PRIVATE ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    target_link_libraries(FloatConservativenessOracleTests PRIVATE interval ${MPFR_LIBRARY} ${GMP_LIBRARY})
    add_test(NAME float_inclusion_mpfr COMMAND FloatConservativenessOracleTests)
    set_tests_properties(float_inclusion_mpfr PROPERTIES LABELS "float-inclusion;mpfr")
    add_executable(FloatBoundsOracleTests tests/float_bounds_tests.cpp)
    target_compile_definitions(FloatBoundsOracleTests PRIVATE INTERVAL_FLOAT_MPFR_ORACLE=1)
    target_include_directories(FloatBoundsOracleTests PRIVATE ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    target_link_libraries(FloatBoundsOracleTests PRIVATE interval ${MPFR_LIBRARY} ${GMP_LIBRARY})
    add_test(NAME float_bounds_mpfr COMMAND FloatBoundsOracleTests)
    set_tests_properties(float_bounds_mpfr PROPERTIES LABELS "float-inclusion;mpfr")

    add_executable(InvalidityOracleTests tests/invalidity_tests.cpp)
    set_target_properties(InvalidityOracleTests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
    target_compile_definitions(InvalidityOracleTests PRIVATE INTERVAL_INVALIDITY_MPFR_ORACLE=1)
    target_include_directories(InvalidityOracleTests PRIVATE ${MPFR_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    target_link_libraries(InvalidityOracleTests PRIVATE interval ${MPFR_LIBRARY} ${GMP_LIBRARY})
    add_test(NAME invalidity_mpfr COMMAND InvalidityOracleTests)
    set_tests_properties(invalidity_mpfr PROPERTIES LABELS "mpfr")
endif()
