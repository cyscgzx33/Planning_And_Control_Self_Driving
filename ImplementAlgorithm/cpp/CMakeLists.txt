project(path_follow)
set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 3.5)

find_package(Qt5Widgets REQUIRED)
add_definitions(-std=c++17)
set(CXX_FLAGS "-Wall -fPIC")
set(CXX_LINKER_FLAGS "-L -Wl,-rpath=/home/honda/Qt/5.12.6/gcc_64/lib -lQt5Core -lQt5Gui")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CXX_LINKER_FLAGS}")

# it's necessary to set RPATH
set(CMAKE_INSTALL_RPATH "/home/honda/Qt/5.12.6/gcc_64/lib")

set(sources path_following_frenet.cpp)

include_directories(${Qt5Widgets_INCLUDES})
add_executable(path_follow ${sources})

#################################### Qt5 library ##################################
add_definitions(${Qt5Widgets_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
#target_link_libraries(path_follow ${Qt5Widgets_LIBRARIES})
###################################################################################

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(path_follow ipopt python2.7 ${Qt5Widgets_LIBRARIES})
else()
target_link_libraries(path_follow ipopt python)
endif()