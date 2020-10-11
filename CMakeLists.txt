cmake_minimum_required(VERSION 3.18.0)
project(input_generator VERSION 0.0.0)

include(CTest)
enable_testing()

set(CMAKE_BUILD_DIRECTORY build)

set (CMAKE_CXX_STANDARD 11)
file(GLOB Voropp_SRC_DIRS "libs/voro++-0.4.6/src")
file(GLOB Blocks_SRC_DIRS "libs/blocks/src")
file(GLOB Boost_LIB_DIRS "/lnm/lib/Q1_2015/boost/lib")

link_directories(${Boost_LIB_DIRS} ${Voropp_SRC_DIRS} ${Blocks_SRC_DIRS})

file(GLOB SOURCES "src/*.cpp")
add_executable(voronoi ${SOURCES})

add_custom_target(
   voropp
   COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/libs/voro++-0.4.6/
)

add_library(blocks STATIC ${Blocks_SRC_DIRS})
set_target_properties(blocks PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(voronoi boost_math_c99l boost_filesystem boost_system voro++ blocks)
add_dependencies(voronoi voropp blocks)
 
# configure_file(CTestCustom.ctest CTestCustom.ctest)
# add_test(straight_filaments_b3r_lin2 bash -c "../run_test.sh straight_filaments_b3r_lin2.sh straight_filaments_b3r_lin2")
# add_test(straight_filaments_b3r_herm2lin3 bash -c "../run_test.sh straight_filaments_b3r_herm2lin3.sh straight_filaments_b3r_herm2lin3")
# add_test(straight_filaments_b3eb bash -c "../run_test.sh straight_filaments_b3eb.sh straight_filaments_b3eb")
# add_test(straight_filaments_b3wk_tan bash -c "../run_test.sh straight_filaments_b3wk_tan.sh straight_filaments_b3wk_tan")
# add_test(voronoi_part6_b3r_herm2lin bash -c "../run_test.sh voronoi_part6_b3r_herm2lin.sh voronoi_b3r_herm2lin3_part6_seed1*")

# Voronoi BC tests (note: you need to choose different seeds, if the rest of the filename is identical!)
# add_test(voronoi_part100_b3r_herm2lin3_dirich bash -c "../run_test.sh voronoi_part100_b3r_herm2lin3_dirich.sh voronoi_b3r_herm2lin3_part100_seed1*")
# add_test(voronoi_part100_b3r_herm2lin3_periodic bash -c "../run_test.sh voronoi_part100_b3r_herm2lin3_periodic.sh voronoi_b3r_herm2lin3_part100_seed2*")
# add_test(voronoi_part100_b3r_herm2lin3_mixed bash -c "../run_test.sh voronoi_part100_b3r_herm2lin3_mixed.sh voronoi_b3r_herm2lin3_part100_seed3*")
# add_test(voronoi_part100_b3r_herm2lin3_thinning bash -c "../run_test.sh voronoi_part100_b3r_herm2lin3_thinning.sh voronoi_b3r_herm2lin3_part100_seed4*")
# add_test(voronoi_part100_b3r_herm2lin3_radDistr bash -c "../run_test.sh voronoi_part100_b3r_herm2lin3_radDistr.sh voronoi_b3r_herm2lin3_part100_seed5*")

# add_test(DNAdoublehelix_2basepairs_b3r_herm2lin3_2eleperturn  bash -c "../run_test.sh DNAdoublehelix_2basepairs_b3r_herm2lin3_2eleperturn.sh DNAdoublehelix_2basepairs_b3r_herm2lin3_2eleperturn")
# add_test(hexagon_b3r_herm2lin3__lfil1_R1e-05 bash -c "../run_test.sh hexagon_b3r_herm2lin3__lfil1_R1e-05.sh hexagon_b3r_herm2lin3__lfil1_R1e-05")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
