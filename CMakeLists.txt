cmake_minimum_required(VERSION 3.16)
project(fem3d_solver VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(fem3d src/fem3d.cpp)
target_include_directories(fem3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(fem3d PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

add_executable(fem3d_solver app/main.cpp)
target_link_libraries(fem3d_solver PRIVATE fem3d)

include(CTest)
if(BUILD_TESTING)
    add_executable(fem3d_tests tests/test_fem3d.cpp)
    target_link_libraries(fem3d_tests PRIVATE fem3d)
    add_test(NAME fem3d_tests COMMAND fem3d_tests)
endif()
