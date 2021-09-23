cmake_minimum_required(VERSION 3.19)
project(geometric_spanners)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")

include_directories(/usr/include)
include_directories(include)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/7/include)

add_executable(geometric_spanners

        include/algorithms/Degree3.h

        include/algorithms/BCC2012.h
        include/algorithms/BGHP2010.h
        include/algorithms/BGS2005.h
        include/algorithms/BHS2018.h
        include/algorithms/BKPX2015.h
        include/algorithms/BSX2009.h
        include/algorithms/KPT2017.h
        include/algorithms/KPX2010.h
        include/algorithms/KX2012.h
        include/algorithms/LW2004.h

        include/printers/GraphPrinter.h
        include/printers/PgfplotPrinter.h
        include/printers/LatexPrinter.h
        include/printers/TablePrinter.h

        include/tools/DelaunayGraph.h
        include/tools/FloydWarshall.h
        include/tools/Metrics.h
        include/tools/PointGenerators.h
        include/tools/Results.h
        include/tools/TDDelaunay.h
        include/tools/Utilities.h

        include/Experiment.h
        include/Names.h
        include/Scratch.h

        main.cpp
        include/algorithms/Degree3.h)

target_link_libraries(geometric_spanners gmp)
target_link_libraries(geometric_spanners mpfr)
target_link_libraries(geometric_spanners boost_system)
target_link_libraries(geometric_spanners boost_thread)
target_link_libraries(geometric_spanners CGAL)
target_link_libraries(geometric_spanners stdc++fs)
