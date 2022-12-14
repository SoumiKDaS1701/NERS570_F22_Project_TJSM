#pragma once
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

#define ERR(e) {std::cout << "Error: " << nc_strerror(e); exit(2);}
using namespace Foam;

class ncdfGrid {
  private:
    int status, ncid, spatial_dim; // Global NetCDF constants, get set by constructor
  public:
    size_t npartitions, max_npr_comm; // Partitioning info
    int npr_comm; // More partitioning info
    size_t nfa_ib, nfa_ibp1; // Face array size info
    size_t nno_ib; // Node array size info
    size_t ncv_ib, ncv_ibp1; // CV array size info
    int nfa_i, nfa_b, nfa_ba, nfa_bi, nfa_bp; // More face array size
    int *pr_comm_list; // List of processors that share faces with current rank
    double* node_cc;
    std::string zone_name[];
    pointField node_xyz; 
    Map<label> ncdfToFoam;
    cellShapeList cells;
    labelList patchToPhys;
    List<DynamicList<face>> patchFaces;
    //polyMesh& mesh;

    /**  
     * @brief Constructor for the ncdfGrid class. Opens the input ncFile and reads in global attributes.
     * @param ncfilein Path to netCDF input file mesh partition
     */
    ncdfGrid(char ncfilein[]){
       this->loadMetaData(ncfilein);
       this->loadPoints();
       this->loadCells();
       this->loadFaces();
    };
   

    // Full List of class methods we want to have
    // 
    // ncdfGrid( )          --- Constructor is complete
    // loadMetaData( )      --- Complete
    // loadPoints( )        --- Complete
    // loadCells( )         --- Complete
    // loadFaces( )         --- Complete
    // writePolyMesh( )     --- Complete in procedural code
    // 
    // loadFlowSolution( )  --- Needs to be done before report
    // interp( )            --- Seems unlikely to be finished before end of semester
    //
    //
    // Other helper/convenience functions:
    //
    // renumber()
    // resize()
    // findFace()
    // findInternalFace()
   
 
    /** 
     * @brief Read the metadata from the grid.nc file
     * @param ncfilein Path to netCDF input file mesh partition
     */
    void loadMetaData(char ncfilein[]){
      
      // NetCDF access variables declaration
      int spatial_dim_id, nprocs_id, max_npr_comm_id, pr_comm_list_id; 
      int ncv_ib_id, ncv_ibp1_id;
      int nfa_i_id, nfa_ib_id, nfa_ibp1_id, nodesOfFace_size_id;
      int nno_ib_id;
      
      
      // open the input file
      Info << "Open the input file" << endl;
      if ((status = nc_open(ncfilein,NC_NOWRITE,&this->ncid)))
          ERR(status);

      // get info on mesh partitioning (1 partition per proc)
      Info << "Get partition info" << endl;
      if ((status = nc_inq_dimid(ncid,"nprocs",&nprocs_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,nprocs_id,&this->npartitions)))
          ERR(status);
      if ((status = nc_get_att_int(ncid,NC_GLOBAL,"npr_comm",&this->npr_comm))) // npr_comm is the number of processors with which this rank must communicate
          ERR(status); 
      if ((status = nc_inq_dimid(ncid,"max_npr_comm",&max_npr_comm_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,max_npr_comm_id,&this->max_npr_comm))) // max_npr_comm is the maximum number of processors with which any rank must communicate
          ERR(status); 
      this->pr_comm_list = (int*)malloc(this->max_npr_comm*sizeof(int));
      if ((status = nc_inq_varid(ncid,"pr_comm_list",&pr_comm_list_id)))
          ERR(status);
      if ((status = nc_get_var_int(ncid,pr_comm_list_id,this->pr_comm_list)))
          ERR(status);
      
      // get the number of faces in current partition
      Info << "Get nfaces info" << endl;
      if ((status = nc_inq_dimid(ncid,"nfa_ib",&nfa_ib_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,nfa_ib_id,&this->nfa_ib)))
          ERR(status);
      if ((status = nc_inq_dimid(ncid,"nfa_ibp1",&nfa_ibp1_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,nfa_ibp1_id,&this->nfa_ibp1)))
          ERR(status);
      if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_i",&this->nfa_i)))
          ERR(status);
      if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_b",&this->nfa_b)))
          ERR(status);
      if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_ba",&this->nfa_ba)))
          ERR(status);
      if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_bi",&this->nfa_bi)))
          ERR(status);
      if ((status = nc_get_att_int(ncid,NC_GLOBAL,"nfa_bp",&this->nfa_bp)))
          ERR(status);
    
      // get the number of cells in current partition
      Info << "Get ncells" << endl;
      if ((status = nc_inq_dimid(ncid,"ncv_ib",&ncv_ib_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,ncv_ib_id,&this->ncv_ib)))
          ERR(status);
      if ((status = nc_inq_dimid(ncid,"ncv_ibp1",&ncv_ibp1_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,ncv_ibp1_id,&this->ncv_ibp1)))
          ERR(status);
   
      // get the number of nodes in current partition and allocate memory
      Info << "Get nnodes" << endl;
      if ((status = nc_inq_dimid(ncid,"nno_ib",&nno_ib_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,nno_ib_id,&this->nno_ib)))
          ERR(status);
    };


 
    /** 
     * @brief Populate the following class attributes from NetCDF grid data:
     *        node_xyz   ----> OpenFOAM pointField object
     *        ncdfToFoam ----> OpenFOAM Map<label> object
     * @param
     */
    void loadPoints(){
      int node_cc_id; // NetCDF variable access 
      this->node_cc = (double*)malloc(3*this->nno_ib*sizeof(double)); 
      label nVerts;
      nVerts = this->nno_ib;
      this->node_xyz.setSize(nVerts);
      this->ncdfToFoam.resize(2*nVerts);
      Info << "Get node array" << endl;
       
      // read in array of node points
      if ((status = nc_inq_varid(ncid,"node_cc",&node_cc_id)))
          ERR(status);
      if ((status = nc_get_var_double(ncid,node_cc_id,this->node_cc)))
          ERR(status);

      Info << "Reshape node points" << endl;
      // Reshape the node points into a 2D array using a Map object
      for(label i=0; i<nVerts; i++){
        point& pt = this->node_xyz[i];
        pt.x() = this->node_cc[3*i];
        pt.y() = this->node_cc[3*i + 1];
        pt.z() = this->node_cc[3*i + 2];
        this->ncdfToFoam.insert(i+1, i);
      }

    };
    
    /** 
     * @brief Populate the cells class attribute. Creates a cellShapeList openFOAM object from NetCDF grid data:
     * @param
     */
    void loadCells(){
      int nodesOfCV_List_id, nodesOfCV_size_id, nodesOfCV_Pointer_id;
      size_t nodesOfCV_size;
      int *nodesOfCV_List, *nodesOfCV_Pointer; 

      Info << "Read connectivity - Nodes of each CV" << endl;

      // Get size of NodesOfCV array
      if ((status = nc_inq_dimid(ncid,"noocv_s",&nodesOfCV_size_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,nodesOfCV_size_id,&nodesOfCV_size)))
          ERR(status);

      // Allocate storage and get the data
      nodesOfCV_List = (int*)malloc(nodesOfCV_size*sizeof(int));
      if ((status = nc_inq_varid(ncid,"noocv_v",&nodesOfCV_List_id)))
          ERR(status);
      if ((status = nc_get_var_int(ncid,nodesOfCV_List_id,nodesOfCV_List)))
          ERR(status);
      // Get pointers for first node of each CV
      nodesOfCV_Pointer = (int*)malloc(this->ncv_ibp1*sizeof(int));
      if ((status = nc_inq_varid(ncid,"noocv_i",&nodesOfCV_Pointer_id)))
        ERR(status);
      if ((status = nc_get_var_int(ncid,nodesOfCV_Pointer_id,nodesOfCV_Pointer)))
        ERR(status);


      labelList hexPoints(8);
      this->cells.setSize(this->ncv_ib);
      const cellModel& hex = cellModel::ref(cellModel::HEX);

      for(label i=0; i<this->ncv_ib; i++){
        hexPoints[0] = nodesOfCV_List[8*i]-1;
        hexPoints[1] = nodesOfCV_List[8*i+1]-1;
        hexPoints[2] = nodesOfCV_List[8*i+2]-1;
        hexPoints[3] = nodesOfCV_List[8*i+3]-1;
        hexPoints[4] = nodesOfCV_List[8*i+4]-1;
        hexPoints[5] = nodesOfCV_List[8*i+5]-1;
        hexPoints[6] = nodesOfCV_List[8*i+6]-1;
        hexPoints[7] = nodesOfCV_List[8*i+7]-1;
        //this->renumber(this->ncdfToFoam, hexPoints);
        this->cells[i] = cellShape(hex, hexPoints);
      } 
    };

    /** 
     * @brief populate openFOAM labelList objects containing face and boundary patch data from NetCDF grid.nc file
     * @param
     */
    void loadFaces(){
      size_t nodesOfFace_size;
      int nzone_id, zone_name_str_len_id, zone_name_id, nodesOfFace_size_id, facesOfZone_size_id;
      int nodesOfFace_List_id, nodesOfFace_Pointer_id, facesOfZone_List_id;
      int *nodesOfFace_List,  *nodesOfFace_Pointer,   *facesOfZone_List;
      size_t nzone, zone_name_str_len; 
      std::string zone_names_str;

      Info << "Read face info" << endl;
      // Get size of FacesOfCV array
      if ((status = nc_inq_dimid(ncid,"noofa_s",&nodesOfFace_size_id)))
          ERR(status);
      if ((status = nc_inq_dimlen(ncid,nodesOfFace_size_id,&nodesOfFace_size)))
          ERR(status);
      // Allocate storage
      nodesOfFace_List = (int*)malloc(nodesOfFace_size*sizeof(int));
      // Get data
      if ((status = nc_inq_varid(ncid,"noofa_v",&nodesOfFace_List_id)))
          ERR(status);
      if ((status = nc_get_var_int(ncid,nodesOfFace_List_id,nodesOfFace_List)))
          ERR(status);
      // Pointers for first face of each CV
      nodesOfFace_Pointer = (int*)malloc(this->nfa_ibp1*sizeof(int));
      if ((status = nc_inq_varid(ncid,"noofa_i",&nodesOfFace_Pointer_id)))
        ERR(status);
      if ((status = nc_get_var_int(ncid,nodesOfFace_Pointer_id,nodesOfFace_Pointer)))
        ERR(status);

      // Boundary info
      Info << "Read boundary zone information" << endl;
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
      if ((status = nc_inq_varid(ncid,"faozn_iv",&facesOfZone_List_id)))
        ERR(status);
      if ((status = nc_get_var_int(ncid,facesOfZone_List_id,facesOfZone_List)))
        ERR(status);
    
      Info << endl;
      Info << "Parsing zone names" << endl;
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
      this->zone_name[nzone];
      zone_names_str = zone_names;
      face quadPoints(4);
      Map<label> physToPatch; // From physical region to Foam patch


      for (label i=0; i<nzone-1; i++)
      {
         Info << endl;
         Info << " Loop nzone: i = " << i << endl;
         this->zone_name[i] = zone_names_str.substr((i+1)*zone_name_str_len,zone_name_str_len);
         this->zone_name[i].erase(std::remove_if(this->zone_name[i].begin(), this->zone_name[i].end(), ::isspace), this->zone_name[i].end());
         Info << "zone_name = " << this->zone_name[i] << endl;

         // label of the face of each zone. end = nfa_ib
         int zone_start_idx = facesOfZone_List[i+1];
         Info << "zone_start_idx = " << zone_start_idx << endl;
         // number of faces in each zone
         int nface_of_zone = 0;
         if (i<nzone-2)
         {
             nface_of_zone = facesOfZone_List[i+2] - facesOfZone_List[i+1];
         }
         else {nface_of_zone = this->nfa_ib - facesOfZone_List[i+1];}
         Info << "nface_of_zone = " << nface_of_zone << endl;

         // Faces of each zone
         for (int j=0; j<nface_of_zone; j++)
         {
           quadPoints[0] = nodesOfFace_List[4*(j + zone_start_idx - 1)]-1;
           quadPoints[1] = nodesOfFace_List[4*(j + zone_start_idx - 1) + 1]-1;
           quadPoints[2] = nodesOfFace_List[4*(j + zone_start_idx - 1) + 2]-1;
           quadPoints[3] = nodesOfFace_List[4*(j + zone_start_idx - 1) + 3]-1;
	   if (j<3){
		   Info << "Before renumber" << endl; 
		   Info << "quadPoints[0] = " << quadPoints[0] << endl; 
		   Info << "nodesOfFace_List[0] = " << nodesOfFace_List[4*(j+zone_start_idx - 1)] << endl; 
	   }
           //Info << "quadPoints assigned j= " << j << endl;
	   //Info << "ncdfToFoam[quadPoints[0]] = " << this->ncdfToFoam[quadPoints[0]] << endl;
	   //Info << "ncdfToFoam[quadPoints[1]] = " << this->ncdfToFoam[quadPoints[1]] << endl;
           //this->renumber(this->ncdfToFoam, quadPoints);
           

	   if (j < 3){
           Info << "quadPoints renumbered " << endl; 
                           Info << "At i,j = " << i << ", " << j << \
                "\n quadPoints      = " << quadPoints[0] << ", " \
                << quadPoints[1] << ", " \
                << quadPoints[2] << ", " \
                << quadPoints[3] << ", "
                
                "\n nodesOfFace_List = " << nodesOfFace_List[4*(j + zone_start_idx - 1)] << ", " \
                << nodesOfFace_List[4*(j + zone_start_idx - 1)+1] << ", " \
                << nodesOfFace_List[4*(j + zone_start_idx - 1)+2] << ", " \
                << nodesOfFace_List[4*(j + zone_start_idx - 1)+3] << endl;
	   }
	   label patchi = -1;
           if (j==0)
           {
               // New region. Allocate patch for it.
               patchi = this->patchFaces.size();
               Info << "patchFaces.size() = " << this->patchFaces.size() << endl;
               this->patchFaces.setSize(patchi + 1);
               this->patchToPhys.setSize(patchi + 1);


               Info<< "Mapping region "<< this->zone_name[i] << " to Foam patch " << patchi << endl;

               physToPatch.insert(i, patchi);
               this->patchToPhys[patchi] = i;
           }
           else { patchi = i; }
           // Add quad to correct patchFaces.
           this->patchFaces[patchi].append(quadPoints);
         } 
      }
    };



    /** 
     * @brief write the openFOAM polyMesh to disk
     * @param
     */
    void writePolyMesh(){
    };
   
 

    /** 
     * @brief renumber an array to openFOAM Map
     * @param mshToFoam const Map<label>
     * @param labels labelList
     */
    void renumber( const Map<label>& mshToFoam, labelList& labels ){
      //Info << "         Entering renumber func" << endl;
      forAll(labels, labelI){
        labels[labelI] = mshToFoam[labels[labelI]];
      }
      //Info << "         Leaving renumber func" << endl;
    };



    /**
      * @brief Find face in pp which uses all vertices in meshF
      * @param pp primitivePatch object
      * @param meshF labelList object to query 
      */
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
       forAll(pFaces, i){
           label facei = pFaces[i];
           const face& f = pp[facei];
           // Count uses of vertices of meshF for f
           label nMatched = 0;
           forAll(f, fp){
               if (meshF.found(f[fp])){
                   nMatched++;
               }
           }

           if (nMatched == meshF.size()){
               return facei;
           }
       }
       return -1;
    };


    /**
      * @brief Find internal face in pp which uses all vertices in meshF
      * @param pp primitivePatch object
      * @param meshF labelList object to query 
      */
    label findInternalFace(const primitiveMesh& mesh, const labelList& meshF){
      const labelList& pFaces = mesh.pointFaces()[meshF[0]];
      forAll(pFaces, i){
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
    };

  
};
