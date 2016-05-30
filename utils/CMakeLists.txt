set(MODULE utils)

set(SRCS
    Utilities.cpp
   )

string(REPLACE ".cpp" ".h" HDRS "${SRCS}")

add_library(${MODULE} ${SRCS})

