/*
 *  geneTree.h
 *  
 *
 *  Created by Jianhong Ou on 7/6/12.
 *  Copyright 2012 UMASSMED. All rights reserved.
 *
 */

#ifndef GENETREE_H_H
#define GENETREE_H_H

#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <deque>
#include <algorithm>
using namespace std;

//node structure;
class node;
//GeneTree class;
class GTree;
//struct for GTree;
struct cmp_ch{
	bool operator()(const char* a, const char* b) const;
};
//enum regulation type;
enum rtype {
	REPRESS=0,//down regulate
	ACTIVE=1,//up regulate
	OSCILLATE=2//others
};

class node{
public:
	const char* name;
	int tol;
	bool miRNA;
	double logFC;
	double pval;
	vector<node*> *par;
	vector<node*> *chd;
public:
	node(const char *t, double l=0.0, bool r=false, double p=0.0);
	~node();
};

//tree
class GTree{
public:
	node *root;
	int size;
	map<const char *, node*, cmp_ch>nodelist;
	bool miRNAcnt4Tol;
	double lfc;
	double pval;
private:
	bool find(const char* s,vector<node*> vec);
	//delete all
	void delall();
	//check foldchange and p value
	bool checkFC(node* n);
	//check direction
	bool checkDir(node* n1, node* n2, rtype d);

public:
	GTree(bool m=false, double l=0.0, double p=0.0, int s=0, node *R=NULL);
	~GTree();
	//search by from
	node* Search(const char* s);
	//change tolerance for children
	void ChangeTolerance(node *n);
	//insert a node into the tree
	void Insert(const char* s, node *n, rtype d);
	//remove a node from the tree
	void Remove(node *n);
	//unbuildPath
	vector<node*> Travel();
	//filter the nodes by lfc
	void verifyFilter(int tolerance = 0);
};

#endif
