#include <gsl/gsl_matrix.h>
#include <list>
#include <vector>
#include <string>
#include "node.h"
#ifndef _tgif_
#define _tgif_
using namespace std;

class TGIF {
	public:
		TGIF(int k, int nitr, int seed, bool verb, double t, double a, bool deb, bool dense, string outputPrefix) : 
				n_components(k),
				max_iter(nitr),
				random_state(seed),
				verbose(verb),
				tol(t),
				alpha(a),
				debug(deb),
				input_dense(dense),
				output_prefix(outputPrefix) {};

		~TGIF() {
			for (vector<Node*>::iterator itr=tree.begin(); itr!=tree.end(); ++itr) {
				gsl_matrix_free((*itr)->get_V());
				if ((*itr)->is_leaf()) {
					gsl_matrix_free(((Leaf*)(*itr))->get_X());
					gsl_matrix_free(((Leaf*)(*itr))->get_U());
				}
			}
			int l = tree.size();
			for (int i = 0; i < l; i++) {
				delete tree[i];
			}
		};

		int fit();
		//int make_tree(vector<int>&, vector<string>&, vector<string>&, vector<int>&, int);
		int make_tree(vector<int>&, vector<string>&, vector<string>&, vector<int>&, int, string, string);
		//int set_same_laplacian_across_tasks(gsl_matrix*, gsl_matrix*);
		int print_factors();

		vector<Node*>* get_tree(){return &tree;};
		list<double>* get_errors(){return &reconstruction_err_;};

	private:
		int n_components;
		int max_iter;
		int random_state;
		bool verbose;
		bool debug;
		double tol;
		double alpha;
		bool input_dense;
		//double lambda;
		vector<Node*> tree;
		list<double> reconstruction_err_;
		string output_prefix;
};
#endif
