cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)

project(nbrfind VERSION 0.1)

add_subdirectory(src)

add_subdirectory(tests)

enable_testing()
