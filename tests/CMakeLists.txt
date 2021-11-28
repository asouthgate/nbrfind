set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(../src/)

add_executable(test_mem test_mem.cpp)

target_link_libraries(test_mem divsufsort uk2 Threads::Threads)

enable_testing()

add_test(TestMem test_mem ${PROJECT_SOURCE_DIR}/res/query_1.fa ${PROJECT_SOURCE_DIR}/res/elan.latest.consensus.16_nov.fasta.h5000)

