cmake_minimum_required(VERSION 3.19)
project(geometric_spanners)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS " -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")

include_directories(/usr/include)
include_directories(include)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/7/include)

add_executable(geometric_spanners
        include/BCC2012.h
        include/BGHP2010.h
        include/BGS2005.h
        include/BHS2017.h
        include/BKPX2015.h
        include/BSX2009.h
        include/DelaunayGraph.h
        include/FloydWarshall.h
        include/GeometricSpannerPrinter.h
        include/GraphAlgoTV.h
        include/KPT2017.h
        include/KPX2010.h
        include/KX2012.h
        include/LW2004.h
        include/metrics.h
        include/TDDelaunay.h
        include/utilities.h
        main.cpp)

target_link_libraries(geometric_spanners gmp)
target_link_libraries(geometric_spanners mpfr)
target_link_libraries(geometric_spanners boost_system)
target_link_libraries(geometric_spanners boost_thread)
target_link_libraries(geometric_spanners CGAL)
