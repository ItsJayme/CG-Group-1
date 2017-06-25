cmake_minimum_required(VERSION 3.7)
project(RayTracing)

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)
#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        imageWriter.h
        main.cpp
        matrix.h
        mesh.cpp
        mesh.h
        RayTracing.cpp
        raytracing.h
        traqueboule.h
        Vec3D.h
        Vertex.h)

add_executable(RayTracing ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(RayTracing ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

install(TARGETS RayTracing RUNTIME DESTINATION bin)
