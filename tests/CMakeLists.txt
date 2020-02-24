project(tests)

# The test program
set(TEST_EXE_DIR ${CMAKE_SOURCE_DIR}/tests/exes)
set(TEST_SRC_DIR ${CMAKE_SOURCE_DIR}/tests/src)

# Add test executables and specify where to put them
add_executable(inputTest ${TEST_SRC_DIR}/inputTest.cpp)
set_target_properties(inputTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(meshTest ${TEST_SRC_DIR}/meshTest.cpp)
set_target_properties(meshTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(materialsTest ${TEST_SRC_DIR}/materialsTest.cpp)
set_target_properties(materialsTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(sgTransportTest ${TEST_SRC_DIR}/sgTransportTest.cpp)
set_target_properties(sgTransportTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(mgTransportTest ${TEST_SRC_DIR}/mgTransportTest.cpp)
set_target_properties(mgTransportTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(sgQuasidiffusionTest ${TEST_SRC_DIR}/sgQuasidiffusionTest.cpp)
set_target_properties(sgQuasidiffusionTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(mgQuasidiffusionTest ${TEST_SRC_DIR}/mgQuasidiffusionTest.cpp)
set_target_properties(mgQuasidiffusionTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(transportQDCouplingTest ${TEST_SRC_DIR}/transportQDCouplingTest.cpp)
set_target_properties(transportQDCouplingTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(qdSolverTest ${TEST_SRC_DIR}/qdSolverTest.cpp)
set_target_properties(qdSolverTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(scbSolverTest ${TEST_SRC_DIR}/scbSolverTest.cpp)
set_target_properties(scbSolverTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

add_executable(startAngleSolverTest ${TEST_SRC_DIR}/startAngleSolverTest.cpp)
set_target_properties(startAngleSolverTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXE_DIR})

# Add the tests
target_link_libraries(inputTest PRIVATE yaml-cpp)
add_test(input ${TEST_EXE_DIR}/inputTest)

target_link_libraries(meshTest PRIVATE libs yaml-cpp)
add_test(mesh ${TEST_EXE_DIR}/meshTest)

target_link_libraries(materialsTest PRIVATE libs yaml-cpp)
add_test(materials ${TEST_EXE_DIR}/materialsTest)

target_link_libraries(mgTransportTest PRIVATE libs yaml-cpp)
add_test(multi_group_transport ${TEST_EXE_DIR}/mgTransportTest)

target_link_libraries(sgTransportTest PRIVATE libs yaml-cpp)
add_test(single_group_transport ${TEST_EXE_DIR}/sgTransportTest)

target_link_libraries(mgQuasidiffusionTest PRIVATE libs yaml-cpp)
add_test(multi_group_quasidiffusion ${TEST_EXE_DIR}/mgQuasidiffusionTest)

target_link_libraries(sgQuasidiffusionTest PRIVATE libs yaml-cpp)
add_test(single_group_quasidiffusion ${TEST_EXE_DIR}/sgQuasidiffusionTest)

target_link_libraries(transportQDCouplingTest PRIVATE libs yaml-cpp)
add_test(transport_qd_coupling ${TEST_EXE_DIR}/transportQDCouplingTest)

target_link_libraries(qdSolverTest PRIVATE libs yaml-cpp)
add_test(quasidiffusion_solver ${TEST_EXE_DIR}/qdSolverTest)

target_link_libraries(startAngleSolverTest PRIVATE libs yaml-cpp)
add_test(starting_angle_solver ${TEST_EXE_DIR}/startAngleSolverTest)

target_link_libraries(scbSolverTest PRIVATE libs yaml-cpp)
add_test(simple_corner_balance_solver ${TEST_EXE_DIR}/scbSolverTest)
