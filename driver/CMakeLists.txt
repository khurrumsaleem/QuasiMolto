add_library(driver driver.cpp)

# Add compilier flags
SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp")
