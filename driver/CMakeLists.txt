add_library(driver driver.cpp)

# Add compilier flags
SET(CMAKE_CXX_FLAGS "-std=c++11 -Ofast -fopenmp -march=native -pipe")
SET(CMAKE_CXX_FLAGS_DEBUG "-pg")
