add_library(uk2 SHARED uk2.cpp probfunc.cpp seqio.cpp sdata.cpp sparseMEM_src/qsufsort.c sparseMEM_src/sparseSA.hpp sparseMEM_src/sparseSA.cpp)

add_library(uk2wrap SHARED uk2wrap.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(nbrfind sparseMEM_src/qsufsort.c sparseMEM_src/sparseSA.hpp sparseMEM_src/sparseSA.cpp main.cpp)

target_link_libraries(nbrfind PRIVATE uk2 Threads::Threads)

add_executable(test_sa sparseMEM_src/qsufsort.c sparseMEM_src/sparseSA.hpp sparseMEM_src/sparseSA.cpp test_sa.cpp)

target_link_libraries(test_sa divsufsort uk2 Threads::Threads)

