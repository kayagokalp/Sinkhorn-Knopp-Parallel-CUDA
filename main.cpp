//#include "scale.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <math.h>       
#include <string.h>
#include <stdlib.h>

#define DEBUG 0 

void wrapper(int* adj, int* xadj, int* tadj, int* txadj, double* rv, double* cv, int* nov, int* nnz, int siter);

void* read_mtxbin(std::string bin_name, int* &adj, int* &xadj, int* &tadj, int* &txadj, double* &rv, double* &cv, int* &nov, int* &nnz){
  
  const char* fname = bin_name.c_str();
  FILE* bp;
  bp = fopen(fname, "rb");
  
  nov = new int;
  nnz = new int;
  
  fread(nov, sizeof(int), 1, bp);
  fread(nnz, sizeof(int), 1, bp);

  std::cout << "READ-nov: " << nov << "\t*nov: " << *nov << "\t&nov: " << &nov << std::endl;
  std::cout << "READ-nnz: " << nnz << "\t*nnz: " << *nnz << "\t&nnz: " << &nnz << std::endl;

  adj = new int[*nnz];   
  xadj = new int[*nov];  
  tadj = new int[*nnz];  
  txadj = new int[*nov];

  fread(adj, sizeof(int), *nnz, bp);
  fread(xadj, sizeof(int), *nov + 1, bp); 

  fread(tadj, sizeof(int), *nnz, bp);
  fread(txadj, sizeof(int), *nov + 1, bp);

  
  std::cout << "Read the binary file" << std::endl;
  
  rv = (double*)malloc(*nov*sizeof(double));
  cv = (double*)malloc(*nov*sizeof(double));
  
  for(int i = 0; i<*nov; i++)
        {
                rv[i] = 1;
                cv[i] = 1;
        }
  
  if(DEBUG){
  std::cout << "##################" << std::endl;
  std::cout << "Binary Read Report" << std::endl;
  std::cout << "nov: " << *nov << std::endl;
  std::cout << "nnz: " << *nnz << std::endl;
  
  for(int i = 0; i < *nov+1; i++){
    std::cout << "i: " << i << "  xadj[i]: " << xadj[i] << std::endl;
  }

  for(int i = 0; i < *nnz; i++){
    std::cout << "i: " << i << "  adj[i]: " << adj[i] << std::endl;
  }

  std::cout << "Binary Read Report" << std::endl;
  std::cout << "##################" << std::endl;
  }
  
}

void fill_man(int* &adj, int* &xadj, int* &tadj, int* &txadj, double* &rv, double* &cv, int* &nov, int* &nnz)
{
	nov = new int;
	nnz = new int;
	*nov = 6;
	*nnz = 9;

	adj = new int[*nnz];   
	xadj = new int[*nov];  
	tadj = new int[*nnz];  
	txadj = new int[*nov];

	adj[0] = 0;
	adj[1] = 2;
	adj[2] = 3;
	adj[3] = 1;
	adj[4] = 1;
	adj[5] = 4;
	adj[6] = 0;
	adj[7] = 2;
	adj[8] = 4;

	xadj[0] = 0;
	xadj[1] = 3;
	xadj[2] = 4;
	xadj[3] = 4;
	xadj[4] = 6;
	xadj[5] = 9;

	tadj[0] = 0;
	tadj[1] = 4;
	tadj[2] = 1;
	tadj[3] = 3;
	tadj[4] = 0;
	tadj[5] = 4;
	tadj[6] = 0;
	tadj[7] = 3;
	tadj[8] = 4;

	txadj[0] = 0;
	txadj[1] = 2;
	txadj[2] = 4;
	txadj[3] = 6;
	txadj[4] = 7;
	txadj[5] = 9;
	
	  
	rv = (double*)malloc(*nov*sizeof(double));
	cv = (double*)malloc(*nov*sizeof(double));

	for(int i = 0; i<*nov; i++)
	{
		rv[i] = 1;
		cv[i] = 1;
	}
if(DEBUG){
  std::cout << "##################" << std::endl;
  std::cout << "Binary Read Report" << std::endl;
  std::cout << "nov: " << *nov << std::endl;
  std::cout << "nnz: " << *nnz << std::endl;
  
  for(int i = 0; i < *nov; i++){
    std::cout << "i: " << i << "  xadj[i]: " << xadj[i] << std::endl;
  }

  for(int i = 0; i < *nnz; i++){
    std::cout << "i: " << i << "  adj[i]: " << adj[i] << std::endl;
  }

  std::cout << "Binary Read Report" << std::endl;
  std::cout << "##################" << std::endl;
  }
}

int main(int argc, char* argv[]){
  
  std::string fname = argv[1];
  std::cout << "fname: " << fname << std::endl;
  int siter = atoi(argv[2]);
  
  int* adj;
  int* xadj;
  int* tadj;
  int* txadj;
  double* rv;
  double* cv;
  int* nov;
  int* nnz;


  //fill_man(adj,xadj,tadj,txadj,rv,cv,nov,nnz);
  read_mtxbin(fname, adj, xadj, tadj, txadj, rv, cv, nov, nnz);
  wrapper(adj, xadj, tadj, txadj, rv, cv, nov, nnz, siter);

}

