# _________________________________________________________________________________
#
# BSD 3-Clause License
#
# Copyright (c) 2021, bionetgen
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# _________________________________________________________________________________
#


cmake_minimum_required(VERSION 3.8.0)
project(input_generator VERSION 0.0.0)

include(CTest)
enable_testing()

set(CMAKE_BUILD_DIRECTORY build)

set(CMAKE_CXX_STANDARD 17)
file(GLOB Voropp_SRC_DIRS "voro++-0.4.6/src")
find_package(Boost COMPONENTS math_c99l REQUIRED)

link_directories(${Voropp_SRC_DIRS})

file(GLOB SOURCES "src/*.cpp")
add_executable(voronoi ${SOURCES})

add_custom_target(
   voropp
   COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/voro++-0.4.6/
)

target_link_libraries(voronoi Boost::math_c99l voro++)
target_include_directories(voronoi PUBLIC ${Boost_INCLUDE_DIRS})
add_dependencies(voronoi voropp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
