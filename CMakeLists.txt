project(Cutting_Plane_Generators)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall")

set(SOURCES main.cpp ./src/cuts/cutter.cpp ./src/parser/parser.cpp)

add_executable(Cutting_Plane_Generators ${SOURCES})
