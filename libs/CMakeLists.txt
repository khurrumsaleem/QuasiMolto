add_library(libs 
	Materials.cpp
        TransportToQDCoupling.cpp 
	MultiGroupQD.cpp 
	SingleGroupQD.cpp 
	Mesh.cpp
        StartingAngle.cpp
        MultiGroupTransport.cpp
        SingleGroupTransport.cpp
	SimpleCornerBalance.cpp
        QuasidiffusionSolver.cpp
	MMS.cpp
        Material.cpp
        SingleGroupDNP.cpp
        MultiGroupDNP.cpp
        MultiGroupQDToMultiPhysicsQDCoupling.cpp
        GreyGroupQD.cpp
        MultiPhysicsCoupledQD.cpp
        HeatTransfer.cpp
        GreyGroupSolver.cpp
        CollapsedCrossSections.cpp
        WriteData.cpp
        MultilevelCoupling.cpp)

target_link_libraries(libs superlu)

# Add compilier flags
#SET(CMAKE_CXX_FLAGS "-std=c++11 -Ofast -fopenmp -pg -msse2")
SET(CMAKE_CXX_FLAGS "-std=c++11 -Ofast -fopenmp")
