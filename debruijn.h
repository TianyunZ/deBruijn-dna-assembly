#pragma once

#ifndef DEBRUIJN_H_INCLUDED
#define DEBRUIJN_H_INCLUDED
#endif

#include "genome.h"
#include "config.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;

struct ArcNode;
struct DBGEdge;

struct DBGNode : Genome {
	int id;
	int in, out;
	vector<DBGNode*> from, to;
	vector<bool> fromVisited, toVisited;
	bool visited;
	ArcNode* firstarc;  // Pointer to the first ArcNode
	DBGNode() {};
	DBGNode(Genome g) {
		clip.assign(g.clip);
		in = 0;
		out = 0;
		visited = false;
		firstarc = NULL;
	};
};

struct ArcNode {
	ArcNode* next;
	DBGNode* node;
	bool visited;
	int repeat;
	bool removed;
	ArcNode(DBGNode* n) {
		node = n;
		next = NULL;
		visited = false;
		repeat = 0;
		removed = false;
	};
};

struct DBGEdge {
	DBGNode* node;
	int cvg;
	int leadToLoop;
	bool visited;
	bool removed;
	DBGEdge(DBGNode* node = nullptr, int cvg = 0)
		: node(node)
		, cvg(cvg)
	{
		leadToLoop = 0;
		visited = false;
		removed = false;
	}
};

struct DeBruijnGraph {  // adjacency list
  private:
    // vector<DBGNode*> adjlist;
    // map<DBGNode, int> nodes;
    vector<DBGNode*> nodes;
    map<Genome, int> IdTable;
  public:
    DeBruijnGraph(){};
    void CreateNode(vector<Genome>&);
    void CreateGraph(vector<Genome>&);
    void EulerianPath(); // Find and return Eulerian path or cycle (as appropriate)

    vector<DBGNode*> FindFirstNode();
    DBGNode* FindNode(Genome);

    void AddEdge(DBGNode*, DBGNode*);
    void DFSHelper(DBGNode*, vector<char>&, DBGNode*);
    void WalkThroughBubble(DBGNode*, vector<DBGNode*>&);
    void RemoveBubble();
	~DeBruijnGraph() {};
};

//void DeBruijnGraph::AddEdge(DBGEdge*& edges, DBGNode* node)
//{
//	if (edges == nullptr) {
//		edges = new DBGEdge(node, 0);
//	}
//	++edges->cvg;
//}

vector<DBGNode*> DeBruijnGraph::FindFirstNode() {
	vector<DBGNode*> res;
	map<Genome, int>::iterator iter;
	cout << IdTable.begin()->second << endl;
	for (iter = IdTable.begin(); iter != IdTable.end(); iter++) {
		/*if ((nodes[iter->second]->out - nodes[iter->second]->in) > 0)*/
		if (nodes[iter->second]->in == 0)
			res.push_back(nodes[iter->second]);
	}
	return res;
}

//void DeBruijnGraph::DFSHelper(DBGEdge* e) {
//	if (e == NULL) {
//		return;
//	}
//	DBGEdge* edge = new DBGEdge();
//	DBGNode* p = new DBGNode();
//	p = edge->node;
//	//DBGNode* last = new DBGNode();
//	vector<char> singleGenome;
//	vector< vector<char> > tempGenomes(4, vector<char>());
//	//singleGenome.push_back(p->clip[0]);
//	for (int n = 0; n < 4; n++) {
//		tempGenomes[n].push_back(p->clip[0]);
//		edge = p->to[n];
//		DFSHelper(edge);
//		//last = edge->node;
//	}
//}

void DeBruijnGraph::DFSHelper(DBGNode* n, vector<char>& singleGenome, DBGNode* last) {
	if (n == NULL) {
		return;
	}
	int i = 0;

	for (i = 0; i < n->to.size(); i ++) {
		if (!n->toVisited[i]) {
			vector<char> tempGenome;
			singleGenome.push_back(n->clip[0]);
			last = n;
			n->toVisited[i] = true;
			DFSHelper(n->to[i], singleGenome, last);
		}
	}
	if (i == n->to.size()) {
		return;
	}
}

void DeBruijnGraph::EulerianPath() {
	vector<vector<char> > resGenome;
	ArcNode* p;
	DBGNode* t;
	int n = 0;
	//DBGEdge* edge = new DBGEdge();
	DBGNode* last = new DBGNode();
	vector<DBGNode*> head = FindFirstNode();
	for (int i = 0; i < head.size(); i++) {
		vector<char> singleGenome;
		//DFSHelper(head[i], singleGenome, last);
		t = head[i];
		while (t->to.size() > 0) {
			for (n = 0; n < t->to.size(); n++) {
				if (!t->toVisited[n]) {
					singleGenome.push_back(t->clip[0]);
					last = t;
					t->toVisited[n] = true;
					t = t->to[n];
					break;
				}
			}
			if (n >= t->to.size()) {
				break;
			}
		}

		//vector< vector<char> > tempGenomes( 4, vector<char>());
		//singleGenome.push_back(head[i]->clip[0]);
		//t = head[i]->to[0];
		//head[i]->toVisited[0] = true;
		//p = head[i]->firstarc;
		/*while (p != NULL) {
			helper(p);
		}
		helper();*/
		/*DFSHelper(head[i]);
		for (int n = 0; n < 4; n ++) {
			edge = head[i]->to[n];
			tempGenomes[n].push_back(edge->node->clip[0]);
			last = edge->node;
			edge = edge->node->to[];
		}*/
		/*while (p != NULL) {
			while (p != NULL && p->visited) {
				p = p->next;
			}
			if (p == NULL) break;
			singleGenome.push_back(p->node->clip[0]);
			p->visited = true;
			last = p->node;
			p = p->node->firstarc;
		}*/
		//while (t != NULL) {
		//	singleGenome.push_back(t->clip[0]);
		//	t->visited = true;
		//	last = t;
		//	for (int n = 0; n < t->to.size(); n++) {
		//		if (!t->toVisited[n]) {
		//			t = t->to[n];
		//			t->toVisited[n] = true;
		//			break;
		//		}
		//	}
		//}
		for (int m = 1; m < last-> clip.length(); m++) {
			singleGenome.push_back(last->clip[m]);
		}
		resGenome.push_back(singleGenome);
		cout << "> short_contig_" << i << endl;
	}

	auto f = fopen(resultPath, "w");
	if (f != NULL)
	{
		int num = 0;
		// fputs("Open File!\n", f);
		for (int i = 0; i < resGenome.size(); i++) {
			if (resGenome[i].size() < minOutPutLength)
				continue;
			fprintf(f, "> short_contig_%d\n", num++);
			for (int j = 0; j < resGenome[i].size(); j++) {
				fprintf(f, "%c", resGenome[i][j]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

DBGNode* DeBruijnGraph::FindNode(Genome mer) {
	auto iter = IdTable.find(mer);
	if (iter == IdTable.end()) {
		auto newNode = new DBGNode(mer);
		newNode->id = nodes.size();
		IdTable[mer] = newNode->id;
		nodes.push_back(newNode);
		return newNode;
	}
	return nodes[iter->second];
}

void DeBruijnGraph::AddEdge(DBGNode* left, DBGNode* right) {
	//ArcNode* newArc = new ArcNode(right);
	//newArc->next = left->firstarc;
	//left->firstarc = newArc;
	//left->out++;
	//right->in++;
	left -> to.push_back(right);
	right -> from.push_back(left);
	left -> toVisited.push_back(false);
	right -> fromVisited.push_back(false);
	left -> out++;
	right -> in++;
}

void DeBruijnGraph::WalkThroughBubble(DBGNode* u, vector<DBGNode*>& path)
{
	ArcNode* p = u->firstarc;
	while (u->out == 1) {
		int i;
		path.push_back(u);
		ArcNode* p = u->firstarc;
		u = p->next->node;
		if (u->in != 1)
			break;
	}
}

//void DeBruijnGraph::RemoveBubble()
//{
//	int numBubble = 0;
//	for (auto&& u : nodes) {
//		ArcNode* x = u->firstarc;
//		while (1) {
//
//		}
//		ArcNode* y;
//		if (x != NULL)
//			ArcNode* y = x->next;
//		if (x != NULL && y != NULL) {
//			if ((double)x->repeat / y->repeat < 0.7) {
//				vector<DBGNode*> xpath, ypath;
//				xpath.push_back(u);
//				ypath.push_back(u);
//				WalkThroughBubble(x->node, xpath);
//				WalkThroughBubble(y->node, ypath);
//
//				ArcNode* tgt_x = xpath.back()->firstarc;
//				ArcNode* tgt_y = ypath.back()->firstarc;
//
//				if (tgt_x != tgt_y || xpath.size() < k / 2) {
//					printf("Walk x, target at %d len = %d cvg = %d\n", tgt_x->id, (int)xpath.size(), x->repeat);
//					for (auto&& t : xpath) {
//						printf("%c", i2c(t.second));
//					}
//					puts("");
//					printf("Walk y, target at %d len = %d cvg = %d\n", tgt_y->id, (int)ypath.size(), y->repeat);
//					for (auto&& t : ypath) {
//						printf("%c", i2c(t.second));
//					}
//					puts("");
//					continue;
//				}
//
//				++numBubble;
//
//				auto lst = xpath.back();
//				auto tgt = lst->to[xpath.back().second]->node;
//				for (int k = 0; k < 4; ++k) {
//					if (tgt->from[k] != nullptr && tgt->from[k]->node == lst) {
//						tgt->from[k] = nullptr;
//					}
//				}
//				u->to[i]->node->removed = true;
//				u->to[i] = nullptr;
//			}
//		}
//
//		for (int i = 0; i < 4; ++i) {
//			for (int j = 0; j < 4; ++j) {
//				auto x = u->to[i];
//				auto y = u->to[j];
//				if (i != j && x != nullptr && y != nullptr) {
//					if ((double)x->cvg / y->cvg < 0.7) {
//						vector<DBGNode*> xpath, ypath;
//						xpath.push_back(u);
//						ypath.push_back(u);
//						WalkThroughBubble(x->node, xpath);
//						WalkThroughBubble(y->node, ypath);
//
//						auto tgt_x = xpath.back().first->to[xpath.back().second]->node;
//						auto tgt_y = ypath.back().first->to[ypath.back().second]->node;
//
//						if (tgt_x != tgt_y || xpath.size() < cfg.k / 2) {
//							printf("Walk x, target at %d len = %d cvg = %d\n", tgt_x->id, (int)xpath.size(), x->cvg);
//							for (auto&& t : xpath) {
//								printf("%c", i2c(t.second));
//							}
//							puts("");
//							printf("Walk y, target at %d len = %d cvg = %d\n", tgt_y->id, (int)ypath.size(), y->cvg);
//							for (auto&& t : ypath) {
//								printf("%c", i2c(t.second));
//							}
//							puts("");
//							continue;
//						}
//
//						++numBubble;
//
//						auto lst = xpath.back().first;
//						auto tgt = lst->to[xpath.back().second]->node;
//						for (int k = 0; k < 4; ++k) {
//							if (tgt->from[k] != nullptr && tgt->from[k]->node == lst) {
//								tgt->from[k] = nullptr;
//							}
//						}
//						u->to[i]->node->removed = true;
//						u->to[i] = nullptr;
//					}
//				}
//			}
//		}
//	}
//	printf("%d bubbles detected\n", numBubble);
//}

void DeBruijnGraph::CreateGraph(vector<Genome>& genomes) // genomes is kmer
{
	printf("Building DBG... ");
	//fflush(stdout);
	for (int i = 0; i < genomes.size(); i++) {
		DBGNode* left = nullptr, * right = nullptr;
		vector<Genome> KmerSet = genomes[i].kmers();
		for (int j = 0; j < KmerSet.size(); j++) {
			left = FindNode(KmerSet[j].leftKm1mer());
			right = FindNode(KmerSet[j].rightKm1mer());
			//AddEdge(left->to[c2i(right->lastChar())], right);
			//AddEdge(right->from[c2i(left->firstChar())], left);
			AddEdge(left, right);
		}
	}
	printf("done\n");
	//fflush(stdout);

	//int t = 0;
	//for (int i = 0; i < genomes.size(); i++) {

	//	DBGNode* leftKmer = new DBGNode(genomes[i]);
	//	leftKmer->clip = genomes[i].clip.substr(0, k-1);
	//	DBGNode* rightKmer = new DBGNode(genomes[i]);
	//	rightKmer->clip = genomes[i].clip.substr(1);
	//	if (nodes.find(*leftKmer) == nodes.end()) {
	//		nodes[*leftKmer] = t++;
	//	}
	//	if (nodes.find(*rightKmer) == nodes.end()) {
	//		//DBGNode* newNode = new DBGNode(genomes[i]);
	//		//adjlist.push_back(newNode);
	//		nodes[*rightKmer] = t++;
	//	}
	//	ArcNode* right = new ArcNode(rightKmer);
	//	leftKmer->firstarc = right;
	//}
}
