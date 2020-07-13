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
    WriteData(Mesh * myMesh, string myOutputDir = "output/");

    // Variables
    string outputDirectory; 
    
    // Pointers 
    Mesh * mesh; 

    // Functions
    void makeDirectory(string myDirName,\
        bool noTimeLabel = false); 
    void write(string myDirName,\
        string parameterName,\
        Eigen::MatrixXd myData,\
        bool noTimeLabel = false);
    void write(string myDirName,\
        string parameterName,\
        Eigen::VectorXd myData,\
        bool noTimeLabel = false);
    void write(string myDirName,\
        string parameterName,\
        double myData,\
        bool noTimeLabel = false);
    void write(string myDirName,\
        string parameterName,\
        int myData,\
        bool noTimeLabel = false);
    string getOutputPath(string myDirName,\
        bool noTimeLabel = false); 

};

//==============================================================================

#endif
