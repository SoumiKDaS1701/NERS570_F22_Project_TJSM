This package provides data processing and conversion between NetCDF and OpenFOAM data types. It uses a WMake based build system to align with the OpenFOAM build framework. It requires that the NetCDF library be installed to the lib/ directory. Additionally, a local copy of OpenFOAM should be installed with sufficient permissions for the user to be able to write a new executable to the `${FOAM_APPBIN}` directory. Before building our executable, the user will need to source the OpenFOAM bashrc to set the relevant `$FOAM_` environment variables. Only one additional configuration step is required for our utility, which is to export an environment variable `${ncdfToFoamRoot}` to the top level directory of this repo as:

    $ export ncdfToFoamRoot=/path/to/git/clone

The ncdfToFoam utility can be compiled from the top level directory of this repo by executing:
       
    $ wmake src
    
The compilation can be verified by entering the test directory and executing:
   
    $ ctest

Final executable should be called as:

    $ ncdfToFoam mesh000.nc

where the present working directory contains:

  
    |-- mesh000.nc
       
    |-- mesh001.nc 
       
    |-- result
  
        |-- mesh000sol.nc
       
        |-- mesh001sol.nc
   
    |-- constant
  
        |-- polyMesh
       
            |-- points
            
            |-- neighbour
            
            |-- owner
            
            |-- faces

    |-- system
        
        |-- controlDict

        |-- fvSolution
