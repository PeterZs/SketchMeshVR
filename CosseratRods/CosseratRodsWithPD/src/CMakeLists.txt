cmake_minimum_required(VERSION 3.5)

#### NOTES ####
#
# compilation with ubuntu reguires at least libglew-dev, freeglut3-dev, libeigen3-dev
#
# sudo apt-get install libglew-dev freeglut3-dev libeigen3-dev
#
# to build and run from this directory call:
#    mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j && ./pd
#
###############

project(ProjectiveDynamicsRods LANGUAGES CXX VERSION 1.0)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# get project files
file(GLOB_RECURSE SRCFILES "projectivedynamics/*.cpp")
add_executable(pd ${SRCFILES})


# include Eigen
FIND_PACKAGE( Eigen3 REQUIRED )
message ("EIGEN FOUND")
target_include_directories( pd PRIVATE ${EIGEN3_INCLUDE_DIR} )

# include/link GLEW
find_package(GLEW REQUIRED)
target_include_directories(pd PRIVATE ${GLEW_INCLUDE_DIRS})
target_link_libraries(pd ${GLEW_LIBRARIES})

# link GLUT
find_package(GLUT REQUIRED)
target_link_libraries(pd ${GLUT_LIBRARIES})

# link GL and GLU
target_link_libraries(pd GL GLU)
