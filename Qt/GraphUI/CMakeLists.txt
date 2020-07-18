cmake_minimum_required(VERSION 3.5)

project(GraphUI LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Boost COMPONENTS system thread REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(CGAL REQUIRED)
#find_package(GEOMETRICSPANNERS REQUIRED)

add_executable(GraphUI
  main.cpp
  mainwindow.cpp
  mainwindow.h
  #mainwindow.ui
  graphview.cpp
  graphview.h
  graphscene.cpp
  graphscene.h
  GraphAlgorithmEvent.h
  DelaunayGraph.h
  GeometricSpannerPrinter.h
  PlanarSpanner.h
  playbackparameters.cpp
  playbackparameters.h
  PolygonSpanner.h
  SpanningGraph.h
  TransformPolygon.h
)

target_link_libraries(GraphUI PRIVATE Qt5::Widgets Boost::system Boost::thread gmp mpfr CGAL::CGAL)
