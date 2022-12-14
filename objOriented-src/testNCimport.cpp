#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <netcdf.h>

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"

#include "IFstream.H"
#include "cellModel.H"
#include "repatchPolyTopoChanger.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointField.H"
#include <ncdfToFoam_test.hpp>

using namespace Foam;

int main(int argc, char *argv[]) {

    argList::addNote
    (
        "Convert NetCDF to OpenFOAM"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word regionName = polyMesh::defaultRegion;

    if (args.readIfPresent("region", regionName))
    {
        Info<< "Creating polyMesh for region " << regionName << endl;
    }

    //  From procedural code main()
    //  // Storage for points
    //  pointField node_xyz;
    //  Map<label> ncdfToFoam;

    //  // Storage for all cells.
    //  cellShapeList cells;
    //  
    //  // Map from patch to gmsh physical region
    //  labelList patchToPhys;
    //  // Storage for patch faces.
    //  List<DynamicList<face>> patchFaces(0);
    //  // Name per physical region
    //  Map<word> physicalNames;

    //  size_t nzone;
    //  size_t zone_name_str_len;
    //  std::string zone_names_str;
    
    
    // Load a test netcdf grid file
    char ncfile[60] = "2x2x2_TestCube.nc";
//    char ncfile[60] = "../Channel-NetCDFgrid/mesh000.nc";
    Info << "Attempting to instantiate ncdfGrid class object from: " << ncfile << endl;
    ncdfGrid partition3(ncfile);
    Info << "Load Complete" << endl;

    //  std::string zone_name[nzone];
    //  
    //  Info << "zone_name_str_len = " << zone_name_str_len << endl;
    //  for (label i=0; i<nzone; i++)
    //  {   
    //      zone_name[i] = zone_names_str.substr(i*zone_name_str_len,zone_name_str_len);
    //      zone_name[i].erase(std::remove_if(zone_name[i].begin(), zone_name[i].end(), ::isspace), zone_name[i].end());
    //      Info << "zone_name = " << zone_name[i] << endl;
    //  }
    
  return 0;
} 
