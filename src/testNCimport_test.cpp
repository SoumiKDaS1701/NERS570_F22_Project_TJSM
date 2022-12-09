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

using namespace Foam;
#define ERR(e) {std::cout << "Error: " << nc_strerror(e) << endl; exit(2);}

void renumber
(
    const Map<label>& mshToFoam,
    labelList& labels
)
{
    forAll(labels, labelI)
    {
        labels[labelI] = mshToFoam[labels[labelI]];
    }
}

// Read in a NetCDF grid from a .nc file
void load
(
    char ncfilein[],
    pointField& node_xyz,
    Map<label>& ncdfToFoam,
    cellShapeList& cells,
    List<DynamicList<face>>& patchFaces,
    
    // Ideally we dont want all of these size_t attributes to be part of the method call signature...
    size_t& npartitions,
    size_t& nno_ib)
{
   
    int ncid, status; // NetCDF library constants
  
    // Grid Partitioning metadata 
    int spatial_dim_id, nprocs_id; // NetCDF variable access ID
    
    // Face arrays
    int    nfa_ib_id, nfa_ibp1_id, nfa_i_id, nodesOfFace_size_id; // NetCDF variable access metadata
    size_t nfa_ib,    nfa_ibp1,              nodesOfFace_size;    // Size metadata
    int                            nfa_i;                         // Size metadata
    int  nodesOfFace_List_id, nodesOfFace_Pointer_id; // NetCDF variable access metadata
    int *nodesOfFace_List,    *nodesOfFace_Pointer;    // Data storage arrays

    // Cell arrays 
    int    ncv_ib_id, ncv_ibp1_id, facesOfCV_size_id, nodesOfCV_size_id;  // NetCDF variable access metadata
    size_t ncv_ib,    ncv_ibp1,    facesOfCV_size,    nodesOfCV_size;     // Size metadata
    int  facesOfCV_List_id, facesOfCV_Pointer_id, nodesOfCV_List_id, nodesOfCV_Pointer_id; // NetCDF variable access metadata
    int *facesOfCV_List,    *facesOfCV_Pointer,    *nodesOfCV_List,    *nodesOfCV_Pointer;    // Data storage arrays

    // Node arrays
    int nno_ib_id, node_cc_id; // NetCDF variable access ID 
    double* node_cc; // data array
    

    // open the input file
    Info << "Open the input file" << endl;
    if ((status = nc_open(ncfilein,NC_NOWRITE,&ncid)))
        ERR(status);

    Info << "Get nprocs" << endl;  
    // get nprocs of partitioned meshes (1 partition per proc)
    if ((status = nc_inq_dimid(ncid,"nprocs",&nprocs_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,nprocs_id,&npartitions)))
        ERR(status);

    Info << "Get nfaces" << endl;      
    // get the number of faces in current partition
    if ((status = nc_inq_dimid(ncid,"nfa_ib",&nfa_ib_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,nfa_ib_id,&nfa_ib)))
        ERR(status);
    if ((status = nc_inq_dimid(ncid,"nfa_ibp1",&nfa_ibp1_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,nfa_ibp1_id,&nfa_ibp1)))
        ERR(status);

    Info << "Get nfaces_internal" << endl;          
    // get the number of internal faces in current partition
    if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_i",&nfa_i)))
        ERR(status);

    Info << "Get ncells" << endl;              
    // get the number of cells in current partition
    if ((status = nc_inq_dimid(ncid,"ncv_ib",&ncv_ib_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,ncv_ib_id,&ncv_ib)))
        ERR(status);
    if ((status = nc_inq_dimid(ncid,"ncv_ibp1",&ncv_ibp1_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,ncv_ibp1_id,&ncv_ibp1)))
        ERR(status);

    Info << "Get nnodes" << endl;                
    // get the number of nodes in current partition and allocate memory
    if ((status = nc_inq_dimid(ncid,"nno_ib",&nno_ib_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,nno_ib_id,&nno_ib)))
        ERR(status);
    Info << "nno_ibp1" << endl; 
    node_cc = (double*)malloc(3*nno_ib*sizeof(double));
    label nVerts;
    nVerts = nno_ib;    
    Info << "Number of Nodes (nVerts) : " << nVerts << endl;
    node_xyz.setSize(nVerts);
    ncdfToFoam.resize(2*nVerts);

    Info << "Get node array" << endl;
    // read in array of node points
    if ((status = nc_inq_varid(ncid,"node_cc",&node_cc_id)))
        ERR(status);
    if ((status = nc_get_var_double(ncid,node_cc_id,node_cc)))
        ERR(status);

    Info << "Reshape node points" << endl;
    // Reshape the node points into a 2D array using a FOAM vector object
    for(label i=0; i<nVerts; i++){
        point& pt = node_xyz[i];
        pt.x() = node_cc[3*i];
        pt.y() = node_cc[3*i + 1];
        pt.z() = node_cc[3*i + 2];

        ncdfToFoam.insert(i+1, i);
    } 

    // Connectivity info - Faces of each CV
    Info << "Read connectivity - Faces of each CV" << endl;
    // Get size of FacesOfCV array
    if ((status = nc_inq_dimid(ncid,"faocv_s",&facesOfCV_size_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,facesOfCV_size_id,&facesOfCV_size)))
        ERR(status);
    Info << "(facesOfCV_size) : " << facesOfCV_size << endl;
    // Allocate storage
    facesOfCV_List = (int*)malloc(facesOfCV_size*sizeof(int));
    // Get data
    if ((status = nc_inq_varid(ncid,"faocv_v",&facesOfCV_List_id)))
        ERR(status);
    if ((status = nc_get_var_int(ncid,facesOfCV_List_id,facesOfCV_List)))
        ERR(status);
    // Pointers for first face of each CV
    facesOfCV_Pointer = (int*)malloc(ncv_ibp1*sizeof(int));
    if ((status = nc_inq_varid(ncid,"faocv_i",&facesOfCV_Pointer_id)))
      ERR(status);
    if ((status = nc_get_var_int(ncid,facesOfCV_Pointer_id,facesOfCV_Pointer)))
      ERR(status);



    // Connectivity info - Nodes of each CV
    Info << "Read connectivity - Nodes of each CV" << endl;
    // Get size of NodesOfCV array
    if ((status = nc_inq_dimid(ncid,"noocv_s",&nodesOfCV_size_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,nodesOfCV_size_id,&nodesOfCV_size)))
        ERR(status);
    Info << "(nodesOfCV_size) : " << nodesOfCV_size << endl;
    //Allocate storage
    nodesOfCV_List = (int*)malloc(nodesOfCV_size*sizeof(int));
    // Get data
    if ((status = nc_inq_varid(ncid,"noocv_v",&nodesOfCV_List_id)))
        ERR(status);
    if ((status = nc_get_var_int(ncid,nodesOfCV_List_id,nodesOfCV_List))) 
        ERR(status);
    for (int i = 0; i < 10; i++){
        Info << "nodesOfCV_List[" << i << "] = " << nodesOfCV_List[i] << endl;
    }
    // Pointers for first node of each CV
    nodesOfCV_Pointer = (int*)malloc(ncv_ibp1*sizeof(int));
    if ((status = nc_inq_varid(ncid,"noocv_i",&nodesOfCV_Pointer_id)))
      ERR(status);
    if ((status = nc_get_var_int(ncid,nodesOfCV_Pointer_id,nodesOfCV_Pointer)))
      ERR(status);



    // Connectivity info - Nodes of each face
    Info << "Read connectivity - Nodes of each face" << endl;
    // Get size of NodesofFace array
    if ((status = nc_inq_dimid(ncid,"noofa_s",&nodesOfFace_size_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,nodesOfFace_size_id,&nodesOfFace_size)))
        ERR(status);
    Info << "(nodesOfFace_size) : " << nodesOfFace_size << endl;
    //Allocate storage
    nodesOfFace_List = (int*)malloc(nodesOfFace_size*sizeof(int));
    if ((status = nc_inq_varid(ncid,"noofa_v",&nodesOfFace_List_id)))
      ERR(status);
    if ((status = nc_get_var_int(ncid,nodesOfFace_List_id,nodesOfFace_List)))
      ERR(status);
    // Pointers for first node of each face
    nodesOfFace_Pointer = (int*)malloc(nfa_ibp1*sizeof(int));
    if ((status = nc_inq_varid(ncid,"noofa_i",&nodesOfFace_Pointer_id)))
      ERR(status);
    if ((status = nc_get_var_int(ncid,nodesOfFace_Pointer_id,nodesOfFace_Pointer)))
      ERR(status);
    


    labelList hexPoints(8);
    cells.setSize(ncv_ib);
    const cellModel& hex = cellModel::ref(cellModel::HEX);
    

    for(label i=0; i<ncv_ib; i++){
        hexPoints[0] = nodesOfCV_List[8*i];
        hexPoints[1] = nodesOfCV_List[8*i+1];
        hexPoints[2] = nodesOfCV_List[8*i+2];
        hexPoints[3] = nodesOfCV_List[8*i+3];
        hexPoints[4] = nodesOfCV_List[8*i+4];
        hexPoints[5] = nodesOfCV_List[8*i+5];
        hexPoints[6] = nodesOfCV_List[8*i+6];
        hexPoints[7] = nodesOfCV_List[8*i+7];
        renumber(ncdfToFoam, hexPoints);
        if (i < 10) { 
            Info << "At i = " << i << \
            "\n hexPoints      = " << hexPoints[0] << ", " \
            << hexPoints[1] << ", " \
            << hexPoints[2] << ", " \
            << hexPoints[3] << ", " \
            << hexPoints[4] << ", " \
            << hexPoints[5] << ", " \
            << hexPoints[6] << ", " \
            << hexPoints[7] << \
            
            "\n nodesOfCV_List = " << nodesOfCV_List[8*i] << ", " \
            << nodesOfCV_List[8*i+1] << ", " \
            << nodesOfCV_List[8*i+2] << ", " \
            << nodesOfCV_List[8*i+3] << ", " \
            << nodesOfCV_List[8*i+4] << ", " \
            << nodesOfCV_List[8*i+5] << ", " \
            << nodesOfCV_List[8*i+6] << ", " \
            << nodesOfCV_List[8*i+7] << endl;
        }
        
        cells[i] = cellShape(hex, hexPoints);
    } 
    for(label i=0; i<10; i++){
        Info << "At i = " << i << ", hexPoints = " << hexPoints[i] << endl;
    } 
    
    face quadPoints(4);
    // From physical region to Foam patch
    Map<label> physToPatch;
    /*
    for (label i=0; i< nodesOfFace_size; i++)
    {
        quadPoints[0] = nodesOfFace_List[4*i];
        quadPoints[1] = nodesOfFace_List[4*i+1];
        quadPoints[2] = nodesOfFace_List[4*i+2];
        quadPoints[3] = nodesOfFace_List[4*i+3];
        
        renumber(ncdfToFoam, quadPoints);
        
        Map<label>::iterator regFnd = physToPatch.find(regPhys);

        label patchi = -1;
        if (regFnd == physToPatch.end())
        {
            // New region. Allocate patch for it.
            patchi = patchFaces.size();

            patchFaces.setSize(patchi + 1);
            patchToPhys.setSize(patchi + 1);

            Info<< "Mapping region " << regPhys << " to Foam patch "
                << patchi << endl;
            physToPatch.insert(regPhys, patchi);
            patchToPhys[patchi] = regPhys;
        }
        else
        {
            // Existing patch for region
            patchi = regFnd();
        }

        // Add quad to correct patchFaces.
        patchFaces[patchi].append(quadPoints);

    }
    

    forAll(patchFaces, patchi)
    {
        patchFaces[patchi].shrink();
    }*/


    // Writing polyMesh/faces file
    Info << "nodesOfFace_List[0] = " << nodesOfFace_List[0] << endl;
    Info << "nodesOfFace_Pointer[0] = " << nodesOfFace_Pointer[0] << endl;
    Info << "nodesOfFace_Pointer[1] = " << nodesOfFace_Pointer[1] << endl;
    Info << nodesOfFace_size << "(" << endl;
    for(int i=0; i<nfa_ibp1-1; i++){
      int NnodesOfThisFace = nodesOfFace_Pointer[i+1] - nodesOfFace_Pointer[i];
      Info << NnodesOfThisFace << "(";
      for(int j=0; j<NnodesOfThisFace; j++){
       Info << nodesOfFace_List[nodesOfFace_Pointer[i] + j - 1] << " ";
      }
      Info << ")" << endl;
    }
    Info << ")" << endl;

};


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

    // Storage for points
    pointField points;
    Map<label> ncdfToFoam;

    // Storage for all cells.
    cellShapeList cells;
    
    
    // Storage for patch faces.
    List<DynamicList<face>> patchFaces(0);

   
    

    // Load a test netcdf grid file
    char ncfile[60] = "../Channel-NetCDFgrid/mesh000.nc";
    Info << "test"<<endl;
    size_t npartitions, nno_ib;
    load(ncfile, points, ncdfToFoam, cells, patchFaces, npartitions, nno_ib);
    //load(ncfile, points, npartitions, nno_ib);


    // Problem is that the orientation of the patchFaces does not have to
    // be consistent with the outwards orientation of the mesh faces. So
    // we have to construct the mesh in two stages:
    // 1. define mesh with all boundary faces in one patch
    // 2. use the read patchFaces to find the corresponding boundary face
    //    and repatch it.
    
    // Create correct number of patches
    // (but without any faces in it)
    faceListList boundaryFaces(patchFaces.size());
    wordList boundaryPatchNames(boundaryFaces.size());
    /*
    forAll(boundaryPatchNames, patchi)
    {
        label physReg = patchToPhys[patchi];

        Map<word>::const_iterator iter = physicalNames.find(physReg);

        if (iter != physicalNames.end())
        {
            boundaryPatchNames[patchi] = iter();
        }
        else
        {
            boundaryPatchNames[patchi] = word("patch") + name(patchi);
        }
        Info<< "Patch " << patchi << " gets name "
            << boundaryPatchNames[patchi] << endl;
    }
    Info<< endl;*/
    
    wordList boundaryPatchTypes(boundaryFaces.size(), polyPatch::typeName);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = polyPatch::typeName;
    wordList boundaryPatchPhysicalTypes
    (
        boundaryFaces.size(),
        polyPatch::typeName
    );
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        std::move(points),
        cells,
        boundaryFaces,
        boundaryPatchNames,
        boundaryPatchTypes,
        defaultFacesName,
        defaultFacesType,
        boundaryPatchPhysicalTypes
    );
    
    
    runTime.setTime(instant(runTime.constant()), 0);
    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    mesh.write();





    // print test results to screen
    std::cout << "num mesh partitions = " << npartitions << std::endl;  
    std::cout << "num nodes in current partitions = " << nno_ib << std::endl;  
    for (size_t i = 0; i < 9; i++)
    {
      std::cout << i << "th point = ("<< points[i].x() << ", " << points[i].y() <<", " << points[i].z() << ")" << std::endl;
    }
    
    std::cout << "# of cells = ("<< cells.size() << ")" << std::endl;
    std::cout << "# of points in a cell = ("<< cells[0].nPoints() << ")" << std::endl;
    std::cout << "# of faces in a cell = ("<< cells[0].nFaces() << ")" << std::endl;
    std::cout << "# of edges in a cell = ("<< cells[0].nEdges() << ")" << std::endl;

    for (size_t i = 0; i < 9; i++)
    {
      //std::cout << i << "th cell = ("<< cells[i] << ")" << std::endl;
    }
  
  return 0;
} 
