/*
 *  geneTree.cpp
 *  
 *
 *  Created by Jianhong Ou on 7/31/12.
 *  Copyright 2012 UMASSMED. All rights reserved.
 *
 */

#include "geneTree.h"
using namespace std;

//node structure;
node::node(const char *t, double l, bool r, double p){
	name=t;
	miRNA=r;
	logFC=l;
	tol=0;
	pval=p;
	par = new vector<node*>;
	chd = new vector<node*>;
}
node::~node(){
	delete par;
	par = NULL;
	delete chd;
	chd = NULL;
}

bool cmp_ch::operator()(const char* a, const char* b){
	return std::strcmp(a,b) < 0 ;
}

//tree private
bool GTree::find(const char* s,vector<node*> vec){
	vector<node*>::iterator it;
	for (it=vec.begin(); it!=vec.end(); it++) {
		if (strcmp((*it)->name,s)==0) {
			return true;
		}
	}
	return false;
}
//delete all
void GTree::delall(){
	map<const char *, node*, cmp_ch>::iterator it;
	for (it=nodelist.begin();it!=nodelist.end();it++) {
		it->second->~node();
	}
	nodelist.erase(nodelist.begin(),nodelist.end());
	delete root;
	root = NULL;
}
//check foldchange and p value
bool GTree::checkFC(node* n){
	bool l = fabs(n->logFC)>=lfc;
	bool p = n->pval<=pval;
	return (l && p);
}
//check direction
bool GTree::checkDir(node* n1, node* n2, rtype d){
	if (d==OSCILLATE || (n1->logFC == 0.0) ) {
		return true;
	}
	if (d==REPRESS) {
		if (checkFC(n1)) {
			if (n1->logFC * n2->logFC < 0.0)
				return true;
			else return false;
		}else {
			return true;
		}
	} else {
		if (d==ACTIVE) {
			if (checkFC(n1)) {
				if (n1->logFC * n2->logFC > 0.0)
					return true;
				else return false;
			} else {
				return true;
			}
		}
	}
	return false;
}

//tree public
GTree::GTree(bool m, double l, double p, int s, node *R):root(R){
	size=s;
	miRNAcnt4Tol=m;
	lfc=l;
	pval=p;
}
GTree::~GTree(){
	if (root) {
		delall();
	}
}
//search by from
node* GTree::Search(const char* s){
	map<const char*, node*>::iterator it;
	it=nodelist.find(s);
	if (it!=nodelist.end()) {
		return it->second;
	} else {
		return NULL;
	}
}
//change tolerance for children
void GTree::ChangeTolerance(node *n){
	if (n->chd) {
		for (int i=0; i<n->chd->size(); i++) {
			int tolerance = checkFC((*n->chd)[i]) ? 0 : 1;
			if ((*n->chd)[i]->miRNA && !miRNAcnt4Tol) tolerance = 0;
			tolerance = n->tol + tolerance;
			if (tolerance < (*n->chd)[i]->tol) {
				(*n->chd)[i]->tol = tolerance;
				ChangeTolerance((*n->chd)[i]);
			}
		}
	}
}
//insert a node into the tree
void GTree::Insert(const char* s, node *n, rtype d){//s is the parent name
	if (root) {
		node* p=Search(s);
		if (p!=NULL && checkDir(p, n, d)) {
			n->tol = checkFC(n) ? 0 : 1;
			if (n->miRNA && !miRNAcnt4Tol) n->tol = 0;
			node* f=Search(n->name);
			if (f!=NULL) {
				(*f->par).push_back(p);
				int tolerance = p->tol + n->tol;
				if (tolerance < f->tol) {
					f->tol = tolerance;
					ChangeTolerance(f);
				}
				(*p->chd).push_back(f);
			} else {
				(*n->par).push_back(p);
				n->tol = p->tol + n->tol;
				(*p->chd).push_back(n);
				nodelist[n->name] = n;
				size++;
			}
		}
	} else {
		root = n;
		size++;
		nodelist[n->name] = n;
	}
}
//remove a node from the tree
void GTree::Remove(node *n){
	if (n->par){
		if (!(*n->par).empty()) {
			vector<node*>::iterator it;
			for (it=(*n->par).begin(); it!=(*n->par).end(); it++) {
				vector<node*> *chd = (*it)->chd;
				vector<node*>::iterator it1;
				for (it1=(*chd).begin(); it1!=(*chd).end(); it1++) {
					if (strcmp((*it1)->name,n->name)==0){
						(*chd).erase(it1);
						break;
					}
				}
			}
		}
	}
	if (n->chd) {
		if (!(*n->chd).empty()) {
			vector<node*>::iterator it;
			for (it=(*n->chd).begin(); it!=(*n->chd).end(); it++) {
				vector<node*> *par = (*it)->par;
				vector<node*>::iterator it1;
				for (it1=(*par).begin(); it1!=(*par).end(); it1++) {
					if (strcmp((*it1)->name,n->name)==0){
						(*par).erase(it1);
						break;
					}
				}
			}
		}
	}
	if(nodelist.find(n->name)!=nodelist.end()) nodelist.erase(nodelist.find(n->name));
	n->~node();
	delete n;
	n=NULL;
	size--;
}
//unbuildPath
vector<node*> GTree::Travel(){
	vector<node*> t;
	deque<node*> Q;
	node* cur;
	if (root) {
		Q.push_back(root);
	} else {
		return t;
	}
	while (!Q.empty()) {
		cur = Q.front();
		if(!find(cur->name,t))
			t.push_back(cur);
		Q.pop_front();
		if (cur->chd) {
			for (int i=0; i<(*cur->chd).size(); i++) {
				if(!find((*cur->chd)[i]->name,t)) {
					Q.push_back((*cur->chd)[i]);
				}
			}
		}
	}
	return t;
}
//filter the nodes by lfc
void GTree::verifyFilter(int tolerance){
	deque<node*> Q;
	vector<node*> t;
	node* cur;
	if (root) {
		Q.push_back(root);
		t.push_back(root);
	} else {
		return ;
	}
	bool changed=false;
	while (!Q.empty()) {
		cur = Q.front();
		Q.pop_front();
		//check and remove
		//is the node without any child?
		bool isLeaf = false;
		bool remove = false;
		if (cur->chd) {
			if (cur->chd->size()==0) {
				isLeaf=true;
			} else {
				if (cur->chd->size()==1) {
					//is the node has only one child and the child point to itself?
					if(strcmp((*cur->chd)[0]->name,cur->name)==0) isLeaf = true;
					else {//is this node is the extra node for other path?
						int next_step_tolerance = checkFC((*cur->chd)[0]) ? 0 : 1;
						if((*cur->chd)[0]->miRNA && !miRNAcnt4Tol) next_step_tolerance = 0;
						if(cur->tol == tolerance && next_step_tolerance) isLeaf = true;
					}					
				}
			}
		} else {
			isLeaf = true;
		}
		vector<node*>::iterator vit=std::find(t.begin(), t.end(), cur);
		if(isLeaf) {
			// remove the node if the logFC==0
			if (!checkFC(cur)) {
				if (strcmp(cur->name,root->name)!=0) {
					if (vit!=t.end()) t.erase(vit);
					Remove(cur);
					remove = true;
				}
				
			}
		} else {
			// remove the node with tolerance > threhold
			if (cur->tol > tolerance) {
				if (strcmp(cur->name,root->name)!=0) {
					if (vit!=t.end()) t.erase(vit);
					Remove(cur);
					remove = true;
				}
			}
		}
		if (!remove) {
			if (cur->chd) {
				for (int i=0; i<(*cur->chd).size(); i++) {
					if(!find((*cur->chd)[i]->name,t)) {
						Q.push_back((*cur->chd)[i]);
						t.push_back((*cur->chd)[i]);
					}
				}
			}
		} else {
			changed=true;
		}
	}
	if (changed) {
		verifyFilter(tolerance);
	}
}
