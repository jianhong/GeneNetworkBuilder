#include <Rcpp.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "geneTree.h"

using namespace Rcpp;
using namespace std;

RcppExport SEXP filterNodes(SEXP xx_from, SEXP xx_to, SEXP xx_miRNA, SEXP xx_logFC, SEXP xx_pval, SEXP xx_dir, SEXP rows, 
							SEXP rootgene, SEXP rootlogFC, SEXP tol, SEXP minify, SEXP miRNAtol, SEXP lFC, SEXP pVAL)
{
	// input the rows
	// xx.from and xx.to must be characters not factors
	// root must be extracted.
	string root=as<string>(rootgene);
	double rootlFC=as<double>(rootlogFC);
	const int nrow=as<const int>(rows);
	bool mini=as<bool>(minify);
	int tolerance=as<int>(tol);
	bool miRNAcntTol=as<bool>(miRNAtol);
	double lfc=as<double>(lFC);
	double pv=as<double>(pVAL);
	//convert the R data.frame into C array
	CharacterVector from(xx_from);
	CharacterVector to(xx_to);
	//convert miRNA into a bool array
	LogicalVector miRNA(xx_miRNA);
	//convert logFC into a double array
	NumericVector logFC(xx_logFC);
	//convert pval into a double array
	NumericVector pval(xx_pval);
	//convert direction into a int array
	IntegerVector dir(xx_dir);
	//buildPath
	GTree gt(miRNAcntTol, lfc, pv);
	//inset rootgene
	node* n= new node(root.c_str(), rootlFC);
	gt.Insert(root.c_str(), n, OSCILLATE);
	//do loop to insert every single line
	vector<int>ids;
	for (int j=0; j<nrow; j++) {
		ids.push_back(j);
	}
	vector<int>::iterator it;
	while(!ids.empty()){
		it=ids.begin();
		while (it!=ids.end()) {
			node* par=gt.Search(from[*it]);
			if (par!=NULL) {
				node* tmp = new node(to[*it],logFC[*it],miRNA[*it],pval[*it]);
				gt.Insert(from[*it],tmp, (rtype)dir[*it]);
				//delete the id
				it=ids.erase(it);
			} else {
				it++;
			}
		}
	}
	//filter the nodes
	gt.verifyFilter(tolerance);
	//out put the nodes
	vector<node*> t=gt.Travel();
	vector<node*>::iterator nit;
	if (t.size()<=1) {
		return List::create(Named(root)=root);
	}
	List result;
	for (nit=t.begin(); nit!=t.end(); nit++) {
		if ((*nit)->par) {
			//calculate the whether the tolerance is come from parent
			int this_step_tolerance = (fabs((*nit)->logFC)>=lfc && (*nit)->pval<=pv) ? 0 : 1;
			if ((*nit)->miRNA && !miRNAcntTol) this_step_tolerance = 0; 
			this_step_tolerance = (*nit)->tol - this_step_tolerance;
			vector<const char*> parents;
			for (unsigned int i=0; i<(*(*nit)->par).size(); i++) {
				if (mini) {
					if ((*(*nit)->par)[i]->tol == this_step_tolerance)
						parents.push_back((*(*nit)->par)[i]->name);
				} else {
					parents.push_back((*(*nit)->par)[i]->name);
				}
			}
			if (parents.size()>0) {
				CharacterVector Rfrom(parents.size());
				for (unsigned int i=0; i<parents.size(); i++) {
					Rfrom[i] = parents[i];
				}
				bool status = false;
				if (mini) {
					if (fabs((*nit)->logFC)<lfc || (*nit)->pval>pv) {
						if ((*nit)->chd) {
							//check whether there is any child come from this node
							for (unsigned int i=0; i<(*(*nit)->chd).size(); i++) {
								int next_step_tolerance = (fabs((*(*nit)->chd)[i]->logFC)>=lfc && (*(*nit)->chd)[i]->pval<=pv) ? 0 : 1;
								if ((*(*nit)->chd)[i]->miRNA && !miRNAcntTol) next_step_tolerance = 0;
								if ((*(*nit)->chd)[i]->tol == (*nit)->tol + next_step_tolerance) {
									status = true;
									break;
								}
							}
						}
					} else {
						status = true;
					}
				} else {
					status = true;
				}
				if(status) result[(*nit)->name] = Rfrom;
			}
		}
	}
	if (result.size()<1) {
		return List::create(Named(root)=root);
	}
	return result;
}
