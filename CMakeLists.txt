cmake_minimum_required (VERSION 2.8.2)
project (CPD-EXAM CXX C)

find_package(OpenCL REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-std=c++11 -D__CL_ENABLE_EXCEPTIONS")

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")

include_directories(${CPD-EXAM_SOURCE_DIR}/event ${CPD-EXAM_SOURCE_DIR}/utils ${CPD-EXAM_SOURCE_DIR}/routines ${OpenCL_INCLUDE_DIRS})
add_subdirectory(event)
add_subdirectory(utils)
add_executable(vertexer main.cc)
target_link_libraries(vertexer utils event ${OpenCL_LIBRARIES})

