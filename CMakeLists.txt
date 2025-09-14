cmake_minimum_required(VERSION 3.16)

project(DiplomSU VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets PrintSupport Concurrent)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets PrintSupport Concurrent)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/form)

set(PROJECT_SOURCES
        src/main.cpp
        src/widget.cpp include/widget.h
        form/widget.ui
        src/landingcalc.cpp include/landingcalc.h
        src/matrix.cpp include/matrix.h
        src/qcustomplot.cpp include/qcustomplot.h
        src/rungekutt.cpp include/rungekutt.h
        src/dynamicRungekutt.cpp include/dynamicRungekutt.h
        src/graphic.cpp include/graphic.h
        src/dynamiccalc.cpp include/dynamiccalc.h
        src/aerodynamic.cpp include/aerodynamic.h
)
#set(app_icon_resource_windows "${CMAKE_CURRENT_SOURCE_DIR}/resources/res.rc")
   # qt_add_executable(photosurfaceexample main.cpp ${app_icon_resource_windows})
qt_add_executable(DiplomSU
       MANUAL_FINALIZATION
       ${PROJECT_SOURCES}
       resources/res.rc
)
target_link_libraries(DiplomSU PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt6::PrintSupport Qt6::Concurrent)


if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(DiplomSU)
endif()
