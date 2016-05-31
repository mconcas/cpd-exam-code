set(MODULE utils)

set(SRCS
    err_code.cpp
    wtime.cpp
   )

   #string(REPLACE ".cpp" ".h" HDRS "${SRCS}")

   add_library(${MODULE} STATIC ${SRCS})

