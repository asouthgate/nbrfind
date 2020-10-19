add_library(uk2 SHARED uk2.cpp probfunc.cpp seqio.cpp)

add_library(uk2wrap SHARED uk2wrap.cpp)

add_executable(nbrfind main.cpp)

target_link_libraries(nbrfind PRIVATE uk2)
