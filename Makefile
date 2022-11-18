ROOT_PATH=/home/tleasca/NERS570_F22_Project_TJSM
LIB=-L$(ROOT_PATH)/lib
INC=-I$(ROOT_PATH)/include -I$(ROOT_PATH)/netcdf-c-4.7.4/include -I$(ROOT_PATH)/src
TEST=test/testNCimport.cpp
SRC=src/interp.cpp
CC=icc

default: test

test:$(TEST)
	$(CC) -Wall -o runTest.exe $(TEST) $(LIB) $(INC) -lnetcdf

main:$(SRC)
	$(CC) -Wall -o interp.exe $(SRC) $(LIB) $(INC) -lnetcdf

clean:
	rm -f ./runTest.exe
