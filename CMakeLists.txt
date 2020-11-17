cmake_minimum_required(VERSION 3.8.0)
project(input_generator VERSION 0.0.0)

include(CTest)
enable_testing()

set(CMAKE_BUILD_DIRECTORY build)

set (CMAKE_CXX_STANDARD 11)
file(GLOB Voropp_SRC_DIRS "voro++-0.4.6/src")
file(GLOB Boost_LIB_DIRS "/lnm/lib/Q1_2015/boost/lib")

link_directories(${Boost_LIB_DIRS} ${Voropp_SRC_DIRS})

file(GLOB SOURCES "src/*.cpp")
add_executable(voronoi ${SOURCES})

add_custom_target(
   voropp
   COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/voro++-0.4.6/
)

target_link_libraries(voronoi boost_math_c99l boost_filesystem boost_system voro++)
add_dependencies(voronoi voropp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
