
# We need to specifically write out this CTestTestfiles.cmake
# because we use wmake instead of cmake

# To complete testing, we have to copy the folder system/ to here
# The folder system/ contains netcdf files to read in

# Source configure in the home directory and wmake in src before testing

# To do tests, type "ctest" in command line


add_subdirectory(/home/tleasca/NERS570_F22_Project_TJSM/src)
add_subdirectory($FOAM_APPBIN)

# Add all that needs to be tested
add_test(ncdfToFoam ncdfToFoam 2x2x2_TestCube.nc)
add_test(checkMesh checkMesh)
