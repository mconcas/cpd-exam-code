cmake_minimum_required (VERSION 2.8.2)
project (CPD-EXAM CXX C)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")

include_directories(${CPD-EXAM_SOURCE_DIR}/base)
add_subdirectory(base)
add_executable(vertexer main.cc)
target_link_libraries(vertexer base)

