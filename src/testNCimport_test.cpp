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
// Find face in pp which uses all vertices in meshF (in mesh point labels)
label findFace(const primitivePatch& pp, const labelList& meshF)
{
    const Map<label>& meshPointMap = pp.meshPointMap();

    // meshF[0] in pp labels.
    if (!meshPointMap.found(meshF[0]))
    {
        Warning<< "Not using gmsh face " << meshF
            << " since zero vertex is not on boundary of polyMesh" << endl;
        return -1;
    }

    // Find faces using first point
    const labelList& pFaces = pp.pointFaces()[meshPointMap[meshF[0]]];

    // Go through all these faces and check if there is one which uses all of
    // meshF vertices (in any order ;-)
    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = pp[facei];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (meshF.found(f[fp]))
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return facei;
        }
    }

    return -1;
}


// Same but find internal face. Expensive addressing.
label findInternalFace(const primitiveMesh& mesh, const labelList& meshF)
{
    const labelList& pFaces = mesh.pointFaces()[meshF[0]];

    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = mesh.faces()[facei];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (meshF.found(f[fp]))
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return facei;
        }
    }
    return -1;
}



// Read in a NetCDF grid from a .nc file
void load
(
    const Foam::string filein,
    pointField& node_xyz,
    Map<label>& ncdfToFoam,
    cellShapeList& cells,
    labelList& patchToPhys,
    List<DynamicList<face>>& patchFaces,
    
    
    size_t& nzone,
    size_t& zone_name_str_len,
    std::string& zone_names_str,
    
    
    // Ideally we dont want all of these size_t attributes to be part of the method call signature...
    size_t& npartitions,
    size_t& nno_ib)
{
    const char* ncfilein = filein.c_str();
   
    int ncid, status; // NetCDF library constants
  
    // Grid Partitioning metadata 
    int spatial_dim_id, nprocs_id, max_npr_comm_id, pr_comm_list_id; // NetCDF variable access ID
    size_t                         max_npr_comm;                     // Size metadata
    int npr_comm, *pr_comm_list;                                     // Data storage

    // Face arrays
    int    nfa_ib_id, nfa_ibp1_id, nfa_i_id, nodesOfFace_size_id; // NetCDF variable access metadata
    size_t nfa_ib,    nfa_ibp1,              nodesOfFace_size;    // Size metadata
    int                            nfa_i, nfa_b, nfa_ba, nfa_bi, nfa_bp;                         // Size metadata
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
    
    // Boundary zone info
    int    nzone_id, zone_name_str_len_id, zone_name_id, facesOfZone_id;    
    int *facesOfZone_List;

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
    if ((status = nc_get_att_int(ncid,NC_GLOBAL,"npr_comm",&npr_comm))) // npr_comm is the number of processors with which myrank communicates
        ERR(status); 
    if ((status = nc_inq_dimid(ncid,"max_npr_comm",&max_npr_comm_id)))
        ERR(status);
    if ((status = nc_inq_dimlen(ncid,max_npr_comm_id,&max_npr_comm))) // max_npr_comm is the maximum number of processors with which any rank communicates
        ERR(status); 
    pr_comm_list = (int*)malloc(max_npr_comm*sizeof(int));
    if ((status = nc_inq_varid(ncid,"pr_comm_list",&pr_comm_list_id)))
        ERR(status);
    if ((status = nc_get_var_int(ncid,pr_comm_list_id,pr_comm_list)))
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
    if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_b",&nfa_b)))
        ERR(status);
    if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_ba",&nfa_ba)))
        ERR(status);
    if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_bi",&nfa_bi)))
        ERR(status);
    if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_bp",&nfa_bp)))
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
    
    // Boundary info
    if ((status = nc_inq_dimid(ncid,"nzone",&nzone_id)))
      ERR(status);
    if ((status = nc_inq_dimlen(ncid,nzone_id,&nzone)))
      ERR(status);
    if ((status = nc_inq_dimid(ncid,"zone_name_str_len",&zone_name_str_len_id)))
      ERR(status);
    if ((status = nc_inq_dimlen(ncid,zone_name_str_len_id,&zone_name_str_len)))
      ERR(status);
    char* zone_names = (char*)malloc(nzone*zone_name_str_len*sizeof(char));
    if ((status = nc_inq_varid(ncid,"zone_name",&zone_name_id)))
      ERR(status);
    if ((status = nc_get_var_text(ncid,zone_name_id,zone_names)))
      ERR(status);
    facesOfZone_List = (int*)malloc(nzone*sizeof(int));
    if ((status = nc_inq_varid(ncid,"faozn_iv",&facesOfZone_id)))
      ERR(status);
    if ((status = nc_get_var_int(ncid,facesOfZone_id,facesOfZone_List)))
      ERR(status);
    Info << "nzone = " << nzone << endl;
    Info << "zone_name_str_len = " << zone_name_str_len << endl;
    Info << "facesOfZone_id = " << facesOfZone_id << endl;
    
    
    Info << endl << "nfa_i = " << nfa_i << endl;
    Info <<         "nfa_b = " << nfa_b << endl;
    Info << "nfa_ib = " << nfa_ib << endl;
    
    Info << "Interior          1 : nfa_i                                    = 1:" << nfa_i << endl;
    Info << "Boundary          nfa_i+1 : nfa_ib                             = " << nfa_i+1 << ":"<< nfa_ib << endl;
    Info << "Assigned boundary nfa_i+1 : nfa_i + nfa_ba                     = " << nfa_i+1 << ":"<< nfa_i+nfa_ba << endl;
    Info << "Interior boundary nfa_i + nfa_ba + 1 : nfa_i + nfa_ba + nfa_bi = " << nfa_i+nfa_ba+1 << ":"<< nfa_i+nfa_ba+nfa_bi << endl;
    Info << "Periodic boundary nfa_i + nfa_ba + nfa_bi + 1 : nfa_ib         = " << nfa_i+nfa_ba+nfa_bi+1 << ":"<< nfa_i+nfa_ba+nfa_bi+nfa_bp << endl;
    
    
    Info << "nfa_ibp1 = " << nfa_ibp1 << endl;
    Info << "nfa_ba = " << nfa_ba << endl;
    Info << "nfa_bi = " << nfa_bi << endl;
    Info << "nfa_bp = " << nfa_bp << endl;
    
    Info << endl;
    for(int i=0;i<nzone;i++){
        for (int j=0 ; j<zone_name_str_len;j++){
            Info << zone_names[i*zone_name_str_len + j];    
        }
        Info << endl;
    }
    Info << endl;
    for(int i=0;i<nzone;i++){
        for (int j=0 ; j<zone_name_str_len;j++){
            Info << zone_names[i*zone_name_str_len + j];    
        }Info << endl;
        Info << "facesOfZone_List["<< i <<"] = " << facesOfZone_List[i] << endl;
    }
     
    labelList hexPoints(8);
    cells.setSize(ncv_ib);
    const cellModel& hex = cellModel::ref(cellModel::HEX);
    
    int zone_start_idx_list[3] = {nfa_i+1, nfa_i + nfa_ba + 1, nfa_i + nfa_ba + nfa_bi + 1};
    int nface_of_zone_list[3] = {nfa_ba, nfa_bi, nfa_bp};
    
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
        
        cells[i] = cellShape(hex, hexPoints);
    }
    
    face quadPoints(4);
    // From physical region to Foam patch
    Map<label> physToPatch;

    std::string zone_name[nzone];
    
    zone_names_str = zone_names;

    for (label i=0; i<nzone-1; i++)
    {   
        zone_name[i] = zone_names_str.substr((i+1)*zone_name_str_len,zone_name_str_len);
        zone_name[i].erase(std::remove_if(zone_name[i].begin(), zone_name[i].end(), ::isspace), zone_name[i].end());
        
        // label of the face of each zone. end = nfa_ib
        int zone_start_idx = facesOfZone_List[i+1];
        // number of faces in each zone
        int nface_of_zone = 0;
        if (i<nzone-2)
        {
            nface_of_zone = facesOfZone_List[i+2] - facesOfZone_List[i+1];
        }
        else {nface_of_zone = nfa_ib - facesOfZone_List[i+1];}
        
        //int zone_start_idx = zone_start_idx_list[i];
        //int nface_of_zone = nface_of_zone_list[i];
        // Faces of each zone
        for (int j=0; j<nface_of_zone; j++)
        {
            quadPoints[0] = nodesOfFace_List[4*(j + zone_start_idx - 1)];
            quadPoints[1] = nodesOfFace_List[4*(j + zone_start_idx - 1) + 1];
            quadPoints[2] = nodesOfFace_List[4*(j + zone_start_idx - 1) + 2];
            quadPoints[3] = nodesOfFace_List[4*(j + zone_start_idx - 1) + 3];
            renumber(ncdfToFoam, quadPoints);
            
            label patchi = -1;
            if (j==0)
            {
                // New region. Allocate patch for it.
                patchi = patchFaces.size();
                patchFaces.setSize(patchi + 1);
                patchToPhys.setSize(patchi + 1);
                
                
                Info<< "Mapping region "<< zone_name[i] << " to Foam patch " << patchi << endl;
                
                physToPatch.insert(i, patchi);
                patchToPhys[patchi] = i;
            }
            else { patchi = i; }
            // Add quad to correct patchFaces.
            patchFaces[patchi].append(quadPoints);
            
        }
        
        
    }
    forAll(patchFaces, patchi)
    {
        patchFaces[patchi].shrink();
    }
    
    Foam::string solnFile = filein;
    solnFile = solnFile.substr(0,solnFile.size()-3);
    solnFile = word("result/") + solnFile + word("sol.nc");
    const char* ncSolnFile = solnFile.c_str();

    
    int ncidSoln, cv_velocity_id;
    double *cv_vel;
    cv_vel = (double*)malloc(3*ncv_ib*sizeof(double));
    if ((status = nc_open(ncSolnFile,NC_NOWRITE,&ncidSoln)))
        ERR(status);
    if ((status = nc_inq_varid(ncidSoln,"cv_velocity",&cv_velocity_id)))
        ERR(status);
    if ((status = nc_get_var_double(ncidSoln,cv_velocity_id,cv_vel)))
        ERR(status);
    Info << "FoamFile" << endl;
    Info << "{" << endl;
    Info << "    version     2.0;" << endl;
    Info << "    format      ascii;" << endl;
    Info << "    class       volVectorField;" << endl;
    Info << "    location    \"1000\";" << endl;
    Info << "    object      U;" << endl;
    Info << "}" << endl;
    Info << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    Info << endl;
    Info << "dimensions      [0 1 -1 0 0 0 0];" << endl;
    Info << endl;
    Info << endl;
    Info << "internalField   nonuniform List<vector>" << endl;
    Info << endl;
    Info << ncv_ib << "(" << endl;
    for(int i=0;i<ncv_ib;i++){ 
      Info << "(" << cv_vel[3*i] <<" "<< cv_vel[3*i + 1] <<" "<< cv_vel[3*i + 2] << ")" << endl;
    }
    Info << ")" << endl;
    
    
    
    

    /*// Writing polyMesh/faces file
    Info << "nodesOfFace_List[0] = " << nodesOfFace_List[0] << endl;
    Info << "nodesOfFace_Pointer[0] = " << nodesOfFace_Pointer[0] << endl;
    Info << "nodesOfFace_Pointer[1] = " << nodesOfFace_Pointer[1] << endl;
    Info << nodesOfFace_size/sizeof(int) << "(" << endl;
    for(int i=0; i<nfa_ibp1-1; i++){
      int NnodesOfThisFace = nodesOfFace_Pointer[i+1] - nodesOfFace_Pointer[i];
      Info << NnodesOfThisFace << "(";
      for(int j=0; j<NnodesOfThisFace; j++){
       Info << nodesOfFace_List[nodesOfFace_Pointer[i] + j - 1] << " ";
      }
      Info << ")" << endl;
    }
    Info << ")" << endl;*/

};


int main(int argc, char *argv[]) {

    argList::addNote
    (
        "Convert NetCDF to OpenFOAM"
    );
    argList::addArgument(".nc file");
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    word regionName = polyMesh::defaultRegion;

    if (args.readIfPresent("region", regionName))
    {
        Info<< "Creating polyMesh for region " << regionName << endl;
    }

    // Storage for points
    pointField node_xyz;
    Map<label> ncdfToFoam;

    // Storage for all cells.
    cellShapeList cells;
    
    // Map from patch to gmsh physical region
    labelList patchToPhys;
    // Storage for patch faces.
    List<DynamicList<face>> patchFaces(0);
    // Name per physical region
    Map<word> physicalNames;

    size_t nzone;
    size_t zone_name_str_len;
    std::string zone_names_str;
    
    
    // Load a test netcdf grid file
    //char ncfile[60] = "./2x2x2_TestCube.nc";
    //char ncfile[60] = "../Channel-NetCDFgrid/mesh000.nc";
    
    Info << "test"<<endl;
    size_t npartitions, nno_ib;
    //load(ncfile, 
    load(args[1],
         node_xyz, 
         ncdfToFoam, 
         cells, 
         patchToPhys, 
         patchFaces, 
         
         nzone,
         zone_name_str_len,
         zone_names_str,
         
         
         npartitions, 
         nno_ib);
         
         
    Info << "Load Complete" << endl;
    std::string zone_name[nzone];
    
    Info << "zone_name_str_len = " << zone_name_str_len << endl;
    for (label i=0; i<nzone-1; i++)
    {   
        zone_name[i] = zone_names_str.substr((i+1)*zone_name_str_len,zone_name_str_len);
        zone_name[i].erase(std::remove_if(zone_name[i].begin(), zone_name[i].end(), ::isspace), zone_name[i].end());
        Info << "zone_name = " << zone_name[i] << endl;
    }
    

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
        std::move(node_xyz),
        cells,
        boundaryFaces,
        boundaryPatchNames,
        boundaryPatchTypes,
        defaultFacesName,
        defaultFacesType,
        boundaryPatchPhysicalTypes
    ); 
    
    // Remove files now, to ensure all mesh files written are consistent.
    mesh.removeFiles();
    
    repatchPolyTopoChanger repatcher(mesh);

    // Now use the patchFaces to patch up the outside faces of the mesh.

    // Get the patch for all the outside faces (= default patch added as last)
    const polyPatch& pp = mesh.boundaryMesh().last();
            
    // Storage for faceZones.
    List<DynamicList<label>> zoneFaces(patchFaces.size());


    // Go through all the patchFaces and find corresponding face in pp.
    forAll(patchFaces, patchi)
    {
       
        const DynamicList<face>& pFaces = patchFaces[patchi];

        forAll(pFaces, i)
        {
            const face& f = pFaces[i];
            // Find face in pp using all vertices of f.
            label patchFacei = findFace(pp, f);
            

            if (patchFacei != -1)
            {
                label meshFacei = pp.start() + patchFacei;

                repatcher.changePatchID(meshFacei, patchi);
            }
            else
            {
                // Maybe internal face? If so add to faceZone with same index
                // - might be useful.
                label meshFacei = findInternalFace(mesh, f);

                if (meshFacei != -1)
                {
                    zoneFaces[patchi].append(meshFacei);
                }
                else
                {
                    WarningInFunction
                        << "Could not match face " << f
                        << " to any of the interior or exterior faces"
                        << " that share the same 0th point" << endl;
                }
            }
        }
        
    }
    Info<< nl;
    
    
    
    
    
    
    runTime.setTime(instant(runTime.constant()), 0);
    
    repatcher.repatch();
    
    // Remove empty defaultFaces
    label defaultPatchID = mesh.boundaryMesh().findPatchID(defaultFacesName);
    if (mesh.boundaryMesh()[defaultPatchID].size() == 0)
    {
        List<polyPatch*> newPatchPtrList((mesh.boundaryMesh().size() - 1));
        label newPatchi = 0;
        forAll(mesh.boundaryMesh(), patchi)
        {
            if (patchi != defaultPatchID)
            {
                const polyPatch& patch = mesh.boundaryMesh()[patchi];

                newPatchPtrList[newPatchi] = patch.clone
                (
                    mesh.boundaryMesh(),
                    newPatchi,
                    patch.size(),
                    patch.start()
                ).ptr();

                newPatchi++;
            }
        }
        repatcher.changePatches(newPatchPtrList);
    }
    
    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    mesh.write();

    Info<< "End\n" << endl;


/*
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
    }*/
  
  return 0;
} 
