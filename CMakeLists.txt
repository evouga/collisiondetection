cmake_minimum_required(VERSION 3.1)
project(collisiondetection)

file(GLOB LIBFILES src/*.cpp)
add_library(${PROJECT_NAME} STATIC ${LIBFILES})
target_include_directories(${PROJECT_NAME} PRIVATE "include" )
target_include_directories(${PROJECT_NAME} PRIVATE "$ENV{EIGEN3_INCLUDE_DIR}" )
install(TARGETS ${PROJECT_NAME} DESTINATION lib)