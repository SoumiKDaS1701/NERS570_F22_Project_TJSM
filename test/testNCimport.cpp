#pragma once
#include <stdio.h>
#include <iostream>
#include <ncdfGrid.hpp>

int main() {
  // Load a test netcdf grid file
  ncdfGrid inpNCgrid;
  char ncfile[60] = "Channel-NetCDFgrid/mesh000.nc"; // this doesn't type match the method arguments
  inpNCgrid.load(ncfile);
  
  // print test results to screen
  std::cout << "num mesh partitions = " << inpNCgrid.npartitions << std::endl;  
  std::cout << "num nodes in current partitions = " << inpNCgrid.nno_ib << std::endl;  
  for (size_t i = 0; i < inpNCgrid.nno_ib; i++)
  {
    std::cout << inpNCgrid.node_cc[i] <<std::endl;
    std::cout << i << std::endl;
  }
  
  return 0;
} 
