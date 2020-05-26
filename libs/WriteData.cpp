/// File: WriteData.cpp
// Purpose: output data 
// Date: May 26, 2019

#include "Mesh.h"
#include "WriteData.h"

using namespace std; 

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,\
    Eigen::DontAlignCols,",","\n");

//==============================================================================
/// WriteData class object constructor
///
/// @param [in] nZ number of axial cells
/// @param [in] nR number of radial cells
WriteData::WriteData(Mesh * myMesh)
{
  // assign pointers
  mesh = myMesh;
};
//==============================================================================

//==============================================================================
/// Make a directory where QM is running 
///
/// @param [in] dirName relative path of directory to be made 
void WriteData::makeDirectory(string myDirName)
{
 
  // Get time at present state
  string dir = getOutputPath(myDirName); 
 
  // Parse system command
  string command = "mkdir -p " + dir;  
 
  // Issue system command to create directory
  system(command.c_str());

};
//===============================================================================

//===============================================================================
/// Write a Eigen::MatrixXd variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    Eigen::MatrixXd myData)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names, and create directory. 
  makeDirectory(myDirName);
  string dir = getOutputPath(myDirName);
  string fileName = dir + parameterName + ".csv";
   
  // Open output file, write myData to it, and close.
  outputFile.open(fileName);
  outputFile << myData.format(CSVFormat) << endl;
  outputFile.close();
 
};
//==============================================================================

//===============================================================================
/// Write a Eigen::VectorXd variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    Eigen::VectorXd myData)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName);
  string dir = getOutputPath(myDirName);
  string fileName = dir + parameterName + ".csv";
  
  // Open output file, write myData to it, and close.
  outputFile.open(fileName);
  outputFile << myData.format(CSVFormat) << endl;
  outputFile.close();
 
};
//==============================================================================

//===============================================================================
/// Write a double variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    double myData)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName);
  string dir = getOutputPath(myDirName);
  string fileName = dir + parameterName + ".csv";
  
  // Open output file, write myData to it, and close.
  outputFile.open(fileName);
  outputFile << myData << endl;
  outputFile.close();
 
};
//==============================================================================

//===============================================================================
/// Write a integer variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    int myData)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName);
  string dir = getOutputPath(myDirName);
  string fileName = dir + parameterName + ".csv";
  
  // Open output file, write myData to it, and close.
  outputFile.open(fileName);
  outputFile << myData << endl;
  outputFile.close();
 
};
//==============================================================================

//==============================================================================
/// Get path of output directory 
///
/// @param [in] dirName relative path of directory to be made 
string WriteData::getOutputPath(string myDirName)
{
 
  // Get time at present state
  string dir = outputDirectory + to_string(mesh->ts[mesh->state]) + "/"\
               + myDirName;

  return dir; 

};
//===============================================================================


