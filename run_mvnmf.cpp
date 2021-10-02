#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <list>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "modules/io.h"
#include "modules/nmf.h" 
#include "modules/utils.h"
#include "modules/node.h"
#include "modules/leaf.h"
#include "modules/initialization.h"
#include "modules/root.h"
#include "modules/branch.h"
#include "modules/tgif.h"

int main(int argc, char **argv)
{
	struct timeval beginTime;
	gettimeofday(&beginTime,NULL);

	struct rusage bUsage;
	getrusage(RUSAGE_SELF,&bUsage);

	const char* treeFile;
	int nFeatures = -1;
	int nComponents = -1;
	
	string outputPrefix = string("");
	int randomState = 1010;
	bool verbose = true;
	bool debug = false;
	int maxIter = 300;
	double tol = 1;
	double alpha = 10;
	//double lambda = 0;
	//int neighborhoodRadius = 2; 
	bool inputDense = false;
	string graphFile = "";
	string uFactorFile = "";
	string vFactorFile = "";
	string usage = string("usage.txt");

	int c;
	while((c = getopt(argc, argv, "o:r:m:t:a:g:u:v:sdhe")) != -1)
		switch (c) {
			case 'o':
				outputPrefix = string(optarg);
				break;
			case 'r':
				randomState = atoi(optarg);
				break;
			case 'm':
				maxIter = atoi(optarg);
				break;
			case 't':
				tol = atof(optarg);
				break;
			case 'a':
				alpha = atof(optarg);
				break;
			/*case 'l':
				lambda = atof(optarg);
				break;*/
			case 'g':
				graphFile = string(optarg);
				break;
			case 'u':
				uFactorFile = string(optarg);
				break;
			case 'v':
				vFactorFile = string(optarg);
				break;
			case 's':
				verbose = false;
				break;
			case 'd':
				debug = true;
				break;
			case 'e':
				inputDense = true;
				break;
			case 'h':
				io::print_usage(usage);
				return 0;
			case '?':
				io::print_usage(usage);
				return 0;
			default:
				io::print_usage(usage);
				return 0;
		}	

	if ((argc - optind) < 3) {
		io::print_usage(usage);
		return 1;
	} else {
		treeFile = argv[optind];
		nFeatures = atoi(argv[optind+1]);
		nComponents = atoi(argv[optind+2]);
	}

	string treeFileName(treeFile);
	vector<int> parentIds, nSamples;
	vector<string> aliases, fileNames;
	io::read_tree(treeFileName, parentIds, aliases, fileNames, nSamples);

	//assuming symmetric matrix, same graph across all tasks
	/*gsl_matrix* W;
	gsl_matrix* L;
	if (lambda > 0) {
		W = gsl_matrix_calloc(nFeatures, nFeatures); 
		L = gsl_matrix_calloc(nFeatures, nFeatures); 
		if (graphFile == "") {
			utils::get_neighborhood_graph(W, neighborhoodRadius);	
		} else {
			io::read_sparse_matrix(graphFile, W);	
		}
		utils::get_laplacian(W, L);
		gsl_matrix_memcpy(W, L);
		gsl_matrix_scale(W, lambda);
		gsl_vector_view diag = gsl_matrix_diagonal(W);
		gsl_vector_add_constant(&diag.vector, 1);
		utils::get_inverse(W);	
	}*/

	TGIF tgif = TGIF(nComponents,maxIter,randomState,verbose,tol,alpha,debug,inputDense,outputPrefix);
	//TGIF tgif = TGIF(nComponents,maxIter,randomState,verbose,tol,alpha,lambda);
	tgif.make_tree(parentIds, aliases, fileNames, nSamples, nFeatures, uFactorFile, vFactorFile);
	//tgif.make_tree(parentIds, aliases, fileNames, nSamples, nFeatures);
	
	/*if (lambda > 0) {
		tgif.set_same_laplacian_across_tasks(L, W);
	}*/
	tgif.fit();
	tgif.print_factors();

	/*if (lambda > 0) {
		gsl_matrix_free(W);
		gsl_matrix_free(L);
	}*/

	struct timeval endTime;
	gettimeofday(&endTime,NULL);

	struct rusage eUsage;
	getrusage(RUSAGE_SELF,&eUsage);

	unsigned long int bt = beginTime.tv_sec;
	unsigned long int et = endTime.tv_sec;

	cout << "Total time elapsed: " << et - bt << " seconds" << endl;

	unsigned long int bu = bUsage.ru_maxrss;
	unsigned long int eu = eUsage.ru_maxrss;
	
	cout << "Memory usage: " << (eu - bu)/1000 << "MB" << endl;
	
	return 0;
}
