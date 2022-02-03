cmake_minimum_required(VERSION 3.19)
project(spanners)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")

## UNCOMMENT FOR RELEASE
#set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")
#set(CMAKE_BUILD_TYPE Release)

include_directories(/usr/include)
include_directories(include)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/7/include)

add_executable(spanners

        algorithms/Degree3.h

        algorithms/BCC2012.h
        algorithms/BGHP2010.h
        algorithms/BGS2005.h
        algorithms/BHS2018.h
        algorithms/BKPX2015.h
        algorithms/BSX2009.h
        algorithms/KPT2017.h
        algorithms/KPX2010.h
        algorithms/KX2012.h
        algorithms/LW2004.h
        algorithms/BoundedDegreePlaneSpanners.h

        tools/printers/GraphPrinter.h
        tools/printers/PgfplotPrinter.h
        tools/printers/LatexPrinter.h
        tools/printers/TablePrinter.h

        tools/Delaunay.h
        tools/DelaunayL2.h
        tools/DelaunayLinf.h
        tools/DelaunayTD.h
        tools/FloydWarshall.h
        tools/Metrics.h
        tools/PointGenerators.h
        tools/Results.h
        tools/StretchFactorExact.h
        tools/StretchFactorExperimental.h
        tools/Utilities.h
#
        Analysis.h
        Experiment.h
        Scratch.h

        main.cpp )

target_link_libraries(spanners gmp)
target_link_libraries(spanners mpfr)
target_link_libraries(spanners boost_system)
target_link_libraries(spanners boost_thread)
target_link_libraries(spanners stdc++fs)
