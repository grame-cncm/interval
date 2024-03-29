cmake_minimum_required(VERSION 3.10)

# set the project name
project(TestInterval)

set(CMAKE_CXX_STANDARD 20)

if (NOTIDY)
else ()
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=-*,clang-*,cppcoreguidelines-*,-cppcoreguidelines-avoid-magic-numbers,performance-*,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,llvm-*,modernize-*,-modernize-use-trailing-return-type,cert-*")
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions(-D_USE_MATH_DEFINES)

include_directories(. ./interval)

# add the executable
add_executable(TestInterval interval/intervalFloatCast.cpp interval/intervalRint.cpp interval/intervalAcos.cpp interval/intervalOr.cpp interval/intervalAsin.cpp interval/intervalLog.cpp interval/intervalDelay.cpp interval/intervalRemainder.cpp interval/intervalRsh.cpp interval/intervalNot.cpp interval/intervalNe.cpp interval/intervalSinh.cpp interval/intervalLsh.cpp interval/intervalNeg.cpp interval/intervalPow.cpp interval/intervalTanh.cpp interval/intervalAsinh.cpp interval/intervalExp.cpp interval/intervalTan.cpp interval/intervalAtanh.cpp interval/intervalEq.cpp interval/intervalMax.cpp interval/intervalMod.cpp interval/intervalLog10.cpp interval/intervalInv.cpp interval/intervalSqrt.cpp interval/intervalIntCast.cpp interval/intervalMul.cpp interval/intervalCosh.cpp interval/intervalCeil.cpp interval/intervalAcosh.cpp interval/intervalGe.cpp interval/intervalAbs.cpp interval/intervalAnd.cpp interval/intervalGt.cpp interval/intervalDiv.cpp interval/intervalXor.cpp interval/intervalSin.cpp interval/intervalCos.cpp interval/intervalAtan2.cpp interval/intervalLt.cpp interval/intervalMem.cpp interval/intervalMin.cpp interval/intervalAdd.cpp interval/intervalAtan.cpp interval/intervalLe.cpp interval/intervalFloor.cpp interval/intervalSub.cpp interval/intervalSub.cpp interval/intervalLabel.cpp interval/interval_algebra.cpp interval/check.cpp interval/bitwiseOperations.cpp main.cpp)


