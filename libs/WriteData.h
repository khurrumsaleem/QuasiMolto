#ifndef WRITEDATA_H
#define WRITEDATA_H 

class Mesh; // forward declaration

using namespace std;

//==============================================================================
//! Class to write out data  

class WriteData
{
  public:
    
    // Constructor
    WriteData(Mesh * myMesh);

    // Variables
    string outputDirectory = "output/";    

    // Pointers 
    Mesh * mesh; 

    // Functions
    void makeDirectory(string myDirName); 
    void write(string myDirName,\
        string parameterName,\
        Eigen::MatrixXd myData);
    void write(string myDirName,\
        string parameterName,\
        Eigen::VectorXd myData);
    void write(string myDirName,\
        string parameterName,\
        double myData);
    void write(string myDirName,\
        string parameterName,\
        int myData);
    string getOutputPath(string myDirName); 

};

//==============================================================================

#endif
