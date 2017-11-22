cmake_minimum_required(VERSION 3.6)
project(src)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
    Mesh.cpp
    Mesh.h
    SketchMesh.cpp
    SketchMesh.h
    Stroke.cpp
    Stroke.h
    SurfaceSmoothing.cpp
    SurfaceSmoothing.h)

add_executable(src ${SOURCE_FILES})