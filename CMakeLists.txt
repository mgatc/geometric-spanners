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

#        include/algorithms/Degree3.h
#
#        include/algorithms/BCC2012.h
#        include/algorithms/BGHP2010.h
#        include/algorithms/BGS2005.h
#        include/algorithms/BHS2018.h
#        include/algorithms/BKPX2015.h
#        include/algorithms/BSX2009.h
#        include/algorithms/KPT2017.h
#        include/algorithms/KPX2010.h
#        include/algorithms/KX2012.h
#        include/algorithms/LW2004.h
#        include/algorithms/BoundedDegreePlaneSpanners.h
#
#        include/printers/GraphPrinter.h
#        include/printers/PgfplotPrinter.h
#        include/printers/LatexPrinter.h
#        include/printers/TablePrinter.h
#
#        include/tools/Delaunay.h
#        include/tools/DelaunayL2.h
#        include/tools/DelaunayLinf.h
#        include/tools/DelaunayTD.h
#        include/tools/FloydWarshall.h
#        include/tools/Metrics.h
#        include/tools/PointGenerators.h
#        include/tools/Results.h
#        include/tools/StretchFactorExact.h
#        include/tools/StretchFactorExperimental.h
#        include/tools/Utilities.h
#
#        include/Analysis.h
#        include/Experiment.h
#        include/Scratch.h

        main.cpp )

target_link_libraries(spanners gmp)
target_link_libraries(spanners mpfr)
target_link_libraries(spanners boost_system)
target_link_libraries(spanners boost_thread)
target_link_libraries(spanners stdc++fs)
