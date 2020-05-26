#ifndef WRITEDATA_H
#define WRITEDATA_H 


using namespace std; 

class Mesh; // forward declaration

//==============================================================================
//! Class to write out data  

class WriteData
{
  public:
    
    // Constructor
    WriteData(Mesh * myMesh); 
    Mesh * mesh; 
};

//==============================================================================

#endif
