// Main Method

#include "debruijn.h"
#include <iostream>
#include <string>
#pragma warning(disable:4996)

using namespace std;

int main()
{
	vector<Genome> fasta1 = ReadFromFasta(fastaPath1);
	vector<Genome> fasta2 = ReadFromFasta(fastaPath2);

  for (int i = 0; i < fasta2.size(); i++) {
		fasta1.push_back(fasta2[i].reverse().complement());
	}

  cout << "Hello World!\n";
	DeBruijnGraph* G = new DeBruijnGraph();

	G->CreateGraph(fasta1);
	G->EulerianPath();

	return 0;
}

