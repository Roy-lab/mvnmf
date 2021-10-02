#include <gsl/gsl_matrix.h>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "node.h"
#include "nmf.h"
#include "leaf.h"
#include "branch.h"
#include "root.h"
#include "io.h"
#include "utils.h"
#include "tgif.h"

int TGIF::fit() {
	double old_error = 0;
	for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			double frobError, regError;
			(*node)->calculate_objective(frobError, regError);
			old_error += (frobError + regError);	
	}
	double old_slope;

	ofstream ofs;
	if (debug) {
		string fname = output_prefix + "track.log";
		ofs.open(fname.c_str());
		ofs << "Itr\t";
		for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			if ((*node)->get_parent() != NULL) {
				ofs << (*node)->get_alias() << "\t";
				if((*node)->is_leaf()) {
					ofs << (*node)->get_alias() << "-parent\t";
				}
			}
		}
		ofs << "total\tslope" << endl;
	}
	
	for (int n_iter =0; n_iter < max_iter; n_iter++){
		for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			(*node)->update();
		}
		double error = 0;
		if (debug) {
			ofs << n_iter+1 << "\t";
		}
		for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			double frobError, regError;
			(*node)->calculate_objective(frobError, regError);
			error += (frobError + regError);
			if (debug && ((*node)->get_parent() != NULL)) {
				if ((*node)->is_leaf()) {
					ofs << frobError << "\t" << regError << "\t";
				} else {
					ofs << regError << "\t";
				}
			}
		}
		reconstruction_err_.push_back(error);
		double slope = old_error - error;
		if (verbose) {
			cout << "Itr " << n_iter +1 << " error = " << error << ", slope = " << slope << endl;
		}
		if (debug) {
			ofs << error << "\t" << slope << endl;
		}
		if (slope < tol) {
			if (verbose) {
				cout << "Converged at iteration " << n_iter << endl;	
			}
			break;
		} else {
			old_error = error;
			old_slope = slope;
		}
	}
	if (debug) {
		ofs.close();
	}
	return 0;
}

int TGIF::make_tree(vector<int>& parentIds, vector<string>& aliases, vector<string>& inputFiles, vector<int>& n, int m, string uFactorFile, string vFactorFile) {
	int nTasks = inputFiles.size();
	int nNodes = parentIds.size();

	//make first leaf node:
	gsl_matrix* X1 = gsl_matrix_calloc(n[0],m);
	if (input_dense) {
		io::read_dense_matrix(inputFiles[0], X1);
	} else {
		io::read_sparse_matrix(inputFiles[0], X1);
	}
	gsl_matrix* U1 = gsl_matrix_calloc(n_components, n[0]);
	gsl_matrix* V1 = gsl_matrix_calloc(n_components, m);
	if ((uFactorFile != "") && (vFactorFile != "")) {
		io::read_dense_matrix(uFactorFile, U1);
		io::read_dense_matrix(vFactorFile, V1);
	} else {
		NMF nmf = NMF(n_components, nndsvd_init, 200, random_state, false, 1);
		nmf.fit(X1, U1, V1);
	}
	Leaf* leaf1 = new Leaf(X1, U1, V1, n_components, alpha, aliases[0]);
	//Leaf* leaf1 = new Leaf(X1, U1, V1, n_components, alpha, lambda, aliases[0]);
	tree.push_back(leaf1); 

	//make rest of nodes:
	//Us and Vs are initialized to the Us and Vs of first leaf node.
	for (int i = 1; i < nNodes; i++) {
		gsl_matrix* V = gsl_matrix_alloc(n_components, m);
		gsl_matrix_memcpy(V, V1);
		if (i < nTasks) {
			gsl_matrix* X = gsl_matrix_calloc(n[i],m);
			if (input_dense) {
				io::read_dense_matrix(inputFiles[i], X);
			} else {
				io::read_sparse_matrix(inputFiles[i], X);
			}
			gsl_matrix* U = gsl_matrix_alloc(n_components, n[i]);
			gsl_matrix_memcpy(U, U1);
			//Leaf* leaf = new Leaf(X, U, V, n_components, alpha, lambda, aliases[i]);
			Leaf* leaf = new Leaf(X, U, V, n_components, alpha, aliases[i]);
			tree.push_back(leaf);
		} else {
			if (parentIds[i] == -1) {
				Root* root = new Root(V, n_components, alpha, aliases[i]);
				tree.push_back(root);	
			}	else {
				Branch* branch = new Branch(V, n_components, alpha, aliases[i]);
				tree.push_back(branch);
			}
		}		
	}

	//link parents and children:
	for (int i =0; i < nNodes; i++) {
		int parentId = parentIds[i]-1;
		if (parentId > 0) {
			Node* child = tree[i];
			Node* parent = tree[parentId]; 
			child->set_parent(parent);	
			parent->add_child(child);
		}
	}
	return 0;
} 

/*
int TGIF::set_same_laplacian_across_tasks(gsl_matrix* L, gsl_matrix* Linv) {
	// Note Linv = (L+I)^(-1). We recycle this again and again for updating
	// U with graph reguluarization (i.e. lambda > 0). We don't have to re-
	// calculate this since U is updated first, and each latent feature of V
	// (each of k rows) is unit-normed after last iteration, so we can assume
	// same v_norm (hence L+I).
	for (vector<Node*>::iterator itr=tree.begin(); itr!=tree.end(); ++itr) {
		if ((*itr)->is_leaf()) {
			((Leaf*)(*itr))->set_laplacian(L, Linv);
		}
	}
}
*/

int TGIF::print_factors() {
	for (vector<Node*>::iterator itr=tree.begin(); itr!=tree.end(); ++itr) {
		(*itr)->write_factors_to_file(output_prefix);
	}
}
