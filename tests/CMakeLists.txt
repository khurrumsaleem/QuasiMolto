project(tests)

# The test program
set(TEST_EXE_DIR ${CMAKE_SOURCE_DIR}/tests/exes)
set(TEST_SRC_DIR ${CMAKE_SOURCE_DIR}/tests/src)

# Add test executables and specify where to put them
add_executable(unitDumbyTest ${TEST_SRC_DIR}/dumbyTest.cpp)
set_target_properties(unitDumbyTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

# Add the tests
target_link_libraries(unitDumbyTest PRIVATE libs)
add_test(test_dumby ${TEST_EXE_DIR}/unitDumbyTest)
