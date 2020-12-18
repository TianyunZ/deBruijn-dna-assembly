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
	DBGNode* from[4], * to[4];
	bool fromVisited[4], toVisited[4];
	bool fromRemoved[4], toRemoved[4];
	bool visited;
	int repeat;
	bool removed;
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
	//vector<DBGNode*> adjlist;
	//map<DBGNode, int> nodes;
	vector<DBGNode*> nodes;
	map<Genome, int> IdTable;
public:
	DeBruijnGraph() {};
	//void CreateNode(vector<Genome>&);
	void CreateGraph(vector<Genome>&);
	void EulerianPath(); // Find and return Eulerian path or cycle (as appropriate)
	vector<DBGNode*> FindFirstNode();
	DBGNode* FindNode(Genome);
	void AddEdge(DBGNode*, DBGNode*);
	void DFSHelper(DBGNode*, vector<char>&, int);
	//void DFSHelper(DBGNode*);
	void WalkThroughBubble(DBGNode*, vector<DBGNode*>&);
	void RemoveBubble();
	void RemoveLoop(DBGNode*, DBGNode*, int);
	~DeBruijnGraph() {};
};

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


void DeBruijnGraph::DFSHelper(DBGNode* n, vector<char>& singleGenome, int index) {
	if (n == NULL) {
		return;
	}
	int i = 0;

	for (i = 0; i < 4; i++) {
		if (n->toRemoved[i]) {
			continue;
		}
		if (!n->toVisited[i]) {
			vector<char> tempGenome;
			singleGenome.push_back(n->clip[0]);
			n->toVisited[i] = true;
			DFSHelper(n->to[i], singleGenome, index + 1);
		}
		else if (n->out == 0) {
			return;
		}
		//else {
		//	int size = singleGenome.size();
		//	for (int j = index; j < size; j++) {
		//		singleGenome.pop_back();
		//	}
		//	return;
		//}
	}
}


void DeBruijnGraph::RemoveLoop(DBGNode* n, DBGNode* last, int index) {
	if (n == NULL) {
		return;
	}
	int i = 0;
	for (i = 0; i < 4; i++) {
		if (n->toRemoved[i]) {
			continue;
		}
		if (!n->toVisited[i]) {
			n->toVisited[i] = true;
			if (n->out > 1) {
				last = n;
				index = i;
			}
			//cout << "node" << n->id << endl;
			RemoveLoop(n->to[i], last, index);
			//cout << "break" << endl;
		}
		else if (index > 0) {
			//n->toRemoved[i] = true;
			last->toRemoved[index] = true;
			//cout << "test" << endl;
		}
	}
}

void DeBruijnGraph::EulerianPath() {
	printf("Building Eulerian Path... \n");
	vector<vector<char>> resGenome;
	ArcNode* p;
	DBGNode* t;
	int n = 0;
	//DBGEdge* edge = new DBGEdge();
	DBGNode* last = new DBGNode();
	vector<DBGNode*> head = FindFirstNode();
	for (int i = 0; i < head.size(); i++) {
		vector<char> singleGenome;
		last = head[i];
		printf("RemoveLoop > %d \n", i);
		RemoveLoop(head[i], last, 0);
		//cout << "End Remove" << endl;
		DFSHelper(head[i], singleGenome, 1);
		/*t = head[i];
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
		}*/
		if (last != NULL) {
			for (int m = 1; m < last->clip.length(); m++) {
				singleGenome.push_back(last->clip[m]);
			}
			resGenome.push_back(singleGenome);
		}
		//cout << "> short_contig_" << i << endl;
	}

	auto f = fopen(resultPath, "w");
	if (f != NULL)
	{
		int num = 0;
		//fputs("Open File!\n", f);
		for (int i = 0; i < resGenome.size(); i++) {
			if (resGenome[i].size() < minOutPutLength)
				continue;
			fprintf(f, "> short_contig_%d\n", num++);
			for (int j = 0; j < resGenome[i].size(); j++) {
				fprintf(f, "%c", resGenome[i][j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
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
	else {
		nodes[iter->second]->repeat++;
	}
	return nodes[iter->second];
}

void DeBruijnGraph::AddEdge(DBGNode* left, DBGNode* right) {
	int l = c2i(left->clip[0]);
	int r = c2i(right->clip[0]);
	left->to[l] = right;
	right->from[r] = left;
	left->toVisited[l] = false;
	right->fromVisited[r] = false;
	left->toRemoved[l] = false;
	right->fromRemoved[r] = false;
	left->out++;
	right->in++;
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

void DeBruijnGraph::RemoveBubble()
{
	int numBubble = 0;
	for (auto&& u : nodes) {
		ArcNode* x = u->firstarc;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				DBGNode* x = u->to[i];
				DBGNode* y = u->to[j];
				if (x != y && x != NULL && y != NULL) {
					if ((double)x->repeat / y->repeat < 0.7) {
						vector<DBGNode*> xpath, ypath;
						xpath.push_back(u);
						ypath.push_back(u);
						WalkThroughBubble(x, xpath);
						WalkThroughBubble(y, ypath);

						ArcNode* tgt_x = xpath.back()->firstarc;
						ArcNode* tgt_y = ypath.back()->firstarc;

						if (tgt_x != tgt_y || xpath.size() < k / 2) {
							printf("Walk x, target at %d len = %d cvg = %d\n", tgt_x->node->id, (int)xpath.size(), x->repeat);
							for (auto&& t : xpath) {
								printf("%c", t->clip[0]);
							}
							puts("");
							printf("Walk y, target at %d len = %d cvg = %d\n", tgt_y->node->id, (int)ypath.size(), y->repeat);
							for (auto&& t : ypath) {
								printf("%c", t->clip[0]);
							}
							puts("");
							continue;
						}

						++numBubble;

						auto lst = xpath.back();
						auto tgt = lst->firstarc->node;
						for (int p = 0; p < 4; ++p) {
							if (tgt->from[p] != nullptr && tgt->from[p] == lst) {
								tgt->from[p] = nullptr;
							}
						}
						u->toRemoved[i] = true;
						/*u->to[i] = nullptr;*/
					}
				}
			}
		}
	}
	printf("%d bubbles detected\n", numBubble);
}


void DeBruijnGraph::CreateGraph(vector<Genome>& genomes) // genomes is kmer
{
	printf("Building DBG... \n");
	//fflush(stdout);
	for (int i = 0; i < genomes.size(); i++) {
		DBGNode* left = nullptr, * right = nullptr;
		vector<Genome> KmerSet = genomes[i].kmers();
		for (int j = 0; j < KmerSet.size(); j++) {
			left = FindNode(KmerSet[j].leftKm1mer());
			right = FindNode(KmerSet[j].rightKm1mer());
			AddEdge(left, right);
		}
	}
	printf("done\n");
}
