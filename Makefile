ROOT_PATH=/home/tleasca/NERS570_F22_Project_TJSM
LIB=-L$(ROOT_PATH)/lib
INC=-I$(ROOT_PATH)/include -I$(ROOT_PATH)/netcdf-c-4.7.4/include -I$(ROOT_PATH)/src
SRC=test/testNCimport.cpp
CC=icc

default: main

main:$(SRC)
	$(CC) -Wall -o runTest.exe $(SRC) $(LIB) $(INC) -lnetcdf
