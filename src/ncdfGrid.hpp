#pragma once
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <netcdf.h>

#define ERR(e) {std::cout << "Error: " << nc_strerror(e); exit(2);}

class ncdfGrid {
  private:
  public:
      size_t npartitions, nno_ib;
      double* node_cc; // How do we make this public attribute an array if we don't know the size yet?

      // Constructor for ncdfGrid class
      ncdfGrid(){ /* Does this need to do anything? */ }

      // Read in a NetCDF grid from a .nc file
      void load(char ncfilein[]){
        int ncid, status;
        int spatial_dim_id, nprocs_id;
        int nno_ib_id, ncv_ib_id, node_cc_id;

        // open the input file
        if ((status = nc_open(ncfilein,NC_NOWRITE,&ncid)))
            ERR(status);
      
        // get nprocs of partitioned meshes (1 partition per proc)
        if ((status = nc_inq_dimid(ncid,"nprocs",&nprocs_id)))
            ERR(status);
        if ((status = nc_inq_dimlen(ncid,nprocs_id,&npartitions)))
            ERR(status);
      
        // get the number of nodes in current partition and allocate memory
        if ((status = nc_inq_dimid(ncid,"nno_ib",&nno_ib_id)))
            ERR(status);
        if ((status = nc_inq_dimlen(ncid,nno_ib_id,&nno_ib)))
            ERR(status);
        this->node_cc = (double*)malloc(3*nno_ib*sizeof(double));
        
        // read in array of node points
        if ((status = nc_inq_varid(ncid,"node_cc",&node_cc_id)))
            ERR(status);
        if ((status = nc_get_var_double(ncid,node_cc_id,node_cc)))
            ERR(status);
          // std::vector<std::vector<double>> node_cc_vec;
    }
};
