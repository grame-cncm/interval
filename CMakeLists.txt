cmake_minimum_required(VERSION 3.12)

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

# Experimental and numerical analysis executable.
add_executable(TestInterval main.cpp)
target_link_libraries(TestInterval PRIVATE interval)

# Fast deterministic regression tests used by CTest and CI.
add_executable(IntervalTests tests/interval_tests.cpp)
target_link_libraries(IntervalTests PRIVATE interval)

enable_testing()
add_test(NAME interval_unit_tests COMMAND IntervalTests)
