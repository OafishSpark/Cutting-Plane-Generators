cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror") # -Wall")

file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS "src/*/*.cpp")

project(CutsHiGHS)

set(HIGHS_DIR /usr/local/lib/cmake/highs)

find_package(HIGHS REQUIRED)
find_package(Threads REQUIRED)

add_executable(CutsHiGHS main.cpp ${SRC_FILES})
target_include_directories(CutsHiGHS PUBLIC ./include)
target_link_libraries(CutsHiGHS highs::highs)
