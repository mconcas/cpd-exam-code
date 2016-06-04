set(MODULE vertex)

set(SRCS
  VertexCandidate.cpp
)

string(REPLACE ".cpp" ".h" HDRS "${SRCS}")

include_directories(${CPD-EXAM_SOURCE_DIR}/vertex)

add_library(${MODULE} ${SRCS})
