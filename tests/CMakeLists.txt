project(tests)
#include_directories(../TPLs/yaml-cpp/include)

# The test program
set(TEST_EXE_DIR ${CMAKE_SOURCE_DIR}/tests/exes)
set(TEST_SRC_DIR ${CMAKE_SOURCE_DIR}/tests/src)

# Add test executables and specify where to put them
add_executable(inputTest ${TEST_SRC_DIR}/inputTest.cpp)
set_target_properties(inputTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

# Add the tests
#target_link_libraries(inputTest PRIVATE libs)
target_link_libraries(inputTest PRIVATE yaml-cpp)
add_test(test_input ${TEST_EXE_DIR}/inputTest)
