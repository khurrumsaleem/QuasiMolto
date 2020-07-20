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
WriteData::WriteData(Mesh * myMesh, string myInputDir)
{
  // assign pointers
  mesh = myMesh;
  
  // Set output root directory 
  outputDirectory = myInputDir;    
};
//==============================================================================

//==============================================================================
/// Make a directory where QM is running 
///
/// @param [in] dirName relative path of directory to be made 
void WriteData::makeDirectory(string myDirName,bool noTimeLabel)
{
 
  // Get time at present state
  string dir = getOutputPath(myDirName,noTimeLabel); 
 
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
    Eigen::MatrixXd myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names, and create directory. 
  makeDirectory(myDirName,noTimeLabel);
  string dir = getOutputPath(myDirName,noTimeLabel);
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
    Eigen::VectorXd myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName,noTimeLabel);
  string dir = getOutputPath(myDirName,noTimeLabel);
  string fileName = dir + parameterName + ".csv";
  
  // Open output file, write myData to it, and close.
  outputFile.open(fileName);
  outputFile << myData.format(CSVFormat) << endl;
  outputFile.close();
 
};
//==============================================================================

//===============================================================================
/// Write a Eigen::VectorXi variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    Eigen::VectorXi myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName,noTimeLabel);
  string dir = getOutputPath(myDirName,noTimeLabel);
  string fileName = dir + parameterName + ".csv";
  
  // Open output file, write myData to it, and close.
  outputFile.open(fileName);
  outputFile << myData.format(CSVFormat) << endl;
  outputFile.close();
 
};
//==============================================================================

//===============================================================================
/// Write a std::vector variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    vector<double> myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  Eigen::VectorXd convertedData;

  // Read std::vector data into an Eigen::VectorXd 
  convertedData.setZero(myData.size());
  for (int iElem = 0; iElem < convertedData.size(); iElem++)
  {
    convertedData(iElem) = myData[iElem];
  }

  write(myDirName,parameterName,convertedData,noTimeLabel);

};
//==============================================================================

//===============================================================================
/// Write a std::vector variable to a file
///
/// @param [in] myDirName relative path of directory to be made 
/// @param [in] parameterName label for variable being written
/// @param [in] myData data to write
void WriteData::write(string myDirName,\
    string parameterName,\
    vector<int> myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  Eigen::VectorXi convertedData;

  // Read std::vector data into an Eigen::VectorXd 
  convertedData.setZero(myData.size());
  for (int iElem = 0; iElem < convertedData.size(); iElem++)
  {
    convertedData(iElem) = myData[iElem];
  }

  write(myDirName,parameterName,convertedData,noTimeLabel);

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
    double myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName,noTimeLabel);
  string dir = getOutputPath(myDirName,noTimeLabel);
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
    int myData,\
    bool noTimeLabel)
{

  // Declare output stream.
  ofstream outputFile;

  // Form directory and file names. 
  makeDirectory(myDirName,noTimeLabel);
  string dir = getOutputPath(myDirName,noTimeLabel);
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
string WriteData::getOutputPath(string myDirName, bool noTimeLabel)
{

  string dir;
  
  if (noTimeLabel)
  {
    // Get name of directory
    dir = outputDirectory + myDirName;
  }
  else
  {
    // Get name of directory at current time
    dir = outputDirectory + to_string(mesh->ts[mesh->state]) + "/"\
               + myDirName;
  }
  return dir; 

};
//===============================================================================


