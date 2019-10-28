# Specify the minimum version for CMake

cmake_minimum_required(VERSION 2.8)

# Project's name

project(QuasiMolto)
# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")
include_directories(
   TPLs/yaml-cpp/include
)
add_subdirectory(driver)
add_subdirectory(libs)
add_subdirectory(TPLs/yaml-cpp)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(QuasiMolto 
	       ${PROJECT_SOURCE_DIR}/driver/driver.cpp 
	       ${PROJECT_SOURCE_DIR}/libs/Materials.cpp
	       ${PROJECT_SOURCE_DIR}/libs/MultiGroupQD.cpp
	       ${PROJECT_SOURCE_DIR}/libs/SingleGroupQD.cpp
	       ${PROJECT_SOURCE_DIR}/libs/Transport.cpp
	       ${PROJECT_SOURCE_DIR}/libs/StartingAngle.cpp
	       ${PROJECT_SOURCE_DIR}/libs/Mesh.cpp)

target_link_libraries(
    QuasiMolto
    PUBLIC
    ${ARMADILLO_LIBRARIES} 
    PRIVATE
    yaml-cpp  # This is the "target" exposed by 
              # yaml-cpp's CMakeLists.txt file
)
