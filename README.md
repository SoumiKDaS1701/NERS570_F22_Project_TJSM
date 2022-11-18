Data processing and conversion of NetCDF and OpenFOAM data types

Final executable should be called as:

$ ./interp.exe. <netCDF_soln_dir> <OpenFOAM_grid_dir> <Interp_out_dir>

where <netCDF_soln_dir> contains:

    |-- grid
  
       |-- mesh000.nc
       
       |-- mesh001.nc 
       
    |-- result
  
       |-- mesh000sol.nc
       
       |-- mesh001sol.nc
       

and <OpenFOAM_grid_dir> contains:

     |-- constant
  
         |-- polyMesh
       
            |-- points
            
            |-- neighbour
            
            |-- owner
            
            |-- faces
