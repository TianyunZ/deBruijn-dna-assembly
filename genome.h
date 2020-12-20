#include "config.h"
#include<vector>
#include<string>
#pragma warning(disable:4996)
using namespace std;

char comp(char c)
{
	switch (c) {
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
}

int c2i(char c)
{
	switch (c) {
	case 'A':
		return 0;
	case 'T':
		return 1;
	case 'G':
		return 2;
	case 'C':
		return 3;
	}
}

struct Genome {  // Genome fragment
	string clip;

	vector<Genome> kmers();
	Genome reverse();
	Genome complement();
	//void km1mer();
	Genome leftKm1mer();
	Genome rightKm1mer();
	Genome() {};
	Genome(string s) {
		this->clip.assign(s);
	}
	~Genome() {};
};

vector<Genome> Genome::kmers() {
	vector<Genome> res;
	if (clip.length() < k) {
		return res;
	}
	for (unsigned int i = 0; i < clip.length() - k + 1; ++i) {
		res.push_back(Genome(clip.substr(i, k)));
	}
	return res;
}

Genome Genome::leftKm1mer() {
	Genome* res = new Genome(clip.substr(0, clip.length() - 1));
	return *res;
}

Genome Genome::rightKm1mer() {
	Genome* res = new Genome(clip.substr(1));
	return *res;
}

Genome Genome::reverse()
{
	Genome t(this->clip);
	std::reverse(t.clip.begin(), t.clip.end());
	return t;
}
Genome Genome::complement()
{
	int l = this->clip.length();
	char* s = new char[l + 1];
	memcpy(s, this->clip.c_str(), l);
	for (int i = 0; i < l; ++i) {
		s[i] = comp(s[i]);
	}
	s[l] = 0;
	Genome t(s);
	return t;
}


vector<Genome> ReadFromFasta(string path) {
	vector<Genome> data;
	printf("reading from %s\n", path.c_str());
	auto f = fopen(path.c_str(), "r");
	char s[222222];
	while (fscanf(f, "%*s\n%s\n", s) != EOF) {
		// string str(s);
		auto g = Genome(s);
		data.push_back(g);
		fscanf(f, "%*s\n%*s\n", s);
	}
	printf("%d genomes read from %s\n", (int)data.size(), path.c_str());
	return data;
}
