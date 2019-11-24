project(path_follow)
set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 2.8.1)

find_package(Qt5Widgets REQUIRED)
add_definitions(-std=c++17)
set(sources path_following_frenet.cpp)

#################################### Qt5 library ##################################
include_directories(${Qt5Widgets_INCLUDES})
add_definitions(${Qt5Widgets_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
###################################################################################

add_executable(path_follow ${sources})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(path_follow ipopt python2.7 ${Qt5Widgets_LIBRARIES})
else()
target_link_libraries(path_follow ipopt python)
endif()