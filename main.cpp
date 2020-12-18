// dna_assembly_zty.cpp 
//
#include "debruijn.h"
//#include "config.h"
#include <iostream>
#include <string>
#pragma warning(disable:4996)
using namespace std;

int main()
{
	vector<Genome> data = ReadFromFasta(shortPath1);
	vector<Genome> data2 = ReadFromFasta(shortPath2);
	for (int i = 0; i < data2.size(); i++) {
		data.push_back(data2[i].reverse().complement());
	}
	cout << "Hello World! \n";
	DeBruijnGraph* G = new DeBruijnGraph();
	G->CreateGraph(data);
	//G->RemoveBubble();
	G->EulerianPath();
	return 0;
}

