#include <iostream>
#include <string>
#include <pcl/io/ply_io.h>

// Function declarations
void usage(const std::string& badarg, const std::string& progName);
std::string getFileExtension(const std::string& fileName);

int main(int argc, char** argv)
{
  // Parse arguments
  std::string input_file;
  std::string output_file_pts;
  std::string output_file_tri;

  // Read arguments and create output files
  if (argc < 2) {
    usage(" input file name not specified.", argv[0]);
  } else {
    cerr << "Using " << argv[1] << " as output file name" << endl;
    input_file = argv[1];
    std::string output_file = getFileWithoutExtension(input_file);
    output_file_pts = output_file + ".pts";
    output_file_tri = output_file + ".tri";
  }
  
}

// Usage
void usage(const std::string& badarg, const std::string& progname)
{
  if(badarg != "") cerr << "Error parsing argument: " << badarg << endl;

  cerr << "Usage: " << progname << " <input file>\n\n";
  cerr << "  e.g., PLYFileReader bunny.ply " << endl;
  ceer << "        will create two files: bunny.pts and bunny.tri" << endl;
  exit(1);
}

std::string getFileWithoutExtension(const std::string& fileName)
{
  if (fileName.find_last_of(".") != std::string::npos)
    return fileName.substr(0, fileName.find_last_of("."));
  return "";
}
