cmake_minimum_required(VERSION 3.19)
project(spanners)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")

## UNCOMMENT FOR RELEASE
#set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")
#set(CMAKE_BUILD_TYPE Release)

include_directories(.)
include_directories(/usr/include)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/7/include)

add_executable(spanners


        src/Results.h
        src/Analysis.h
        src/Experiment.h
        src/Scratch.h

        main.cpp)

target_link_libraries(spanners gmp)
target_link_libraries(spanners mpfr)
target_link_libraries(spanners boost_system)
target_link_libraries(spanners boost_thread)
#target_link_libraries(spanners stdc++fs)
