set(MODULE base)

set(SRCS
  Event.cpp
  VertexCandidate.cpp
  VertexerFast.cpp
)

string(REPLACE ".cpp" ".h" HDRS "${SRCS}")


add_library(${MODULE} ${SRCS})
