set(MODULE base)

set(SRCS
  Event.cpp
)

string(REPLACE ".cpp" ".h" HDRS "${SRCS}")

include_directories(${CPD-EXAM_SOURCE_DIR}/utils)

add_library(${MODULE} ${SRCS})
