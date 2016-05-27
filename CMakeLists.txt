cmake_minimum_required (VERSION 2.8.2)
project (CPD-EXAM CXX C)

find_package(OpenMP REQUIRED)
if ( ${CMAKE_C_COMPILER} MATCHES "icc.*$" )
  set(OFFLOAD_FLAGS "${OFFLOAD_FLAGS} -offload-option,mic,ld,\"-rpath /opt/intel/mic/myo/lib/ -rpath /opt/intel/mic/coi/device-linux-release/lib/\"")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${OFFLOAD_FLAGS} -std=c++11")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -opt-report-phase=offloa -opt-report-phase=offload")

include_directories(${CPD-EXAM_SOURCE_DIR}/event ${CPD-EXAM_SOURCE_DIR}/utils ${CPD-EXAM_SOURCE_DIR}/routines)
add_subdirectory(event)
add_executable(vertexer main.cc)
target_link_libraries(vertexer base)

