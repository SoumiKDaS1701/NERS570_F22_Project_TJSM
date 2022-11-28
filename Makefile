CC=icc
ROOT_PATH=/home/tleasca/NERS570_F22_Project_TJSM
OF_ROOT_PATH=$(ROOT_PATH)/OpenFOAM/OpenFOAM-v1812
LIB=-L$(ROOT_PATH)/lib
INC=-I$(ROOT_PATH)/src -I$(ROOT_PATH)/include -I$(ROOT_PATH)/netcdf-c-4.7.4/include -I$(OF_ROOT_PATH)/src/OpenFOAM/lnInclude -I$(OF_ROOT_PATH)/src/OSspecific/POSIX -I$(OF_ROOT_PATH)/src/finiteVolume/fvMesh -I$(OF_ROOT_PATH)/src/finiteVolume/fvMesh/fvBoundaryMesh -I$(OF_ROOT_PATH)/src/finiteVolume/fvMesh/fvPatches/fvPatch/ -I$(OF_ROOT_PATH)/src/finiteVolume/fields/fvPatchFields/fvPatchField/ -I$(OF_ROOT_PATH)/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolation/ -I$(OF_ROOT_PATH)/src/finiteVolume/fields/volFields -I$(OF_ROOT_PATH)/src/finiteVolume/fields/surfaceFields -I$(OF_ROOT_PATH)/src/finiteVolume/finiteVolume/fvSchemes/ -I$(OF_ROOT_PATH)/src/finiteVolume/finiteVolume/fvSolution/


TEST=test/testNCimport.cpp
SRC=src/interp.cpp

default:test

test:$(TEST)
	$(CC) -Wall -g -DWM_DP -DWM_LABEL_SIZE=${WM_LABEL_SIZE} -o runTest.exe $(TEST) $(LIB) $(INC) -lnetcdf
	#$(CC) -Wall -g -o runTest.exe $(TEST) $(LIB) $(INC) -lnetcdf

main:$(SRC)
	$(CC) -Wall -o interp.exe $(SRC) $(LIB) $(INC) -lnetcdf

clean:
	rm -f ./runTest.exe
