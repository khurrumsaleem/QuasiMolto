// File: MultiGroupQDToMultiPhysicsQDCoupling.cpp     
// Purpose: Calculate coupling factors between multi-group quasidiffusion and
//   multi-physics quasidiffusion.
// Date: April 9, 2020

#include "MultiGroupQDToMultiPhysicsQDCoupling.h"

using namespace std;

//==============================================================================
/// MultiGroupQDToMultiPhysicsQDCoupling class object constructor
///
/// @param [in] myMesh mesh object 
/// @param [in] myMats materials object 
/// @param [in] myInput input object 
/// @param [in] myMPQD MultiPhysicsCoupledQD object 
/// @param [in] myMGPQD multigroup quasidiffusion object 
MGQDToMPQDCoupling::MGQDToMPQDCoupling(Mesh * myMesh,\
        Materials * myMats,\
        YAML::Node * myInput,\
        MultiPhysicsCoupledQD * myMPQD,\
        MultiGroupQD * myMGQD)
{

  // Assign inputs to their member variables
  mesh = myMesh; 
  mats = myMats;
  input = myInput;
  mpqd = myMPQD;
  mgqd = myMGQD;

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MGQDToMPQDCoupling::collapseNuclearData()
{

  // Temporary accumulator variables
  double fluxAccum,rCurrentAccum,zCurrentAccum;

  // Reset collapsed nuclear data
  mats->oneGroupXS->resetData();

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;
      rCurrentAccum = 0.0;
      zCurrentAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Flux weighted sigT
        // Flux weighted sigS
        // Flux weighted neutron velocity
        // Flux weighted Ezz
        // Flux weighted Err
        // Flux weighted Erz
        // Flux weighted quasidiffusion coefficient
        // Flux weighted DNP coefficient 
        // Accumulate flux

        // Axial current weighted neutron velocity
        // Axial current weighted sigTR
        // Accumulate axial current

        // Radial current weighted neutron velocity
        // Radial current weighted sigTR
        // Accumulate radial current

        // Radial current weighted zeta1
        // Radial current weighted zeta2

        // Axial current weighted zeta1
        // Axial current weighted zeta2
      }
    }
  }

};
//==============================================================================
