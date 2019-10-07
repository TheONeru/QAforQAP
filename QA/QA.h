#pragma once
#include "common.h"
class QA
{
public:
	QA(parameters p, int **D, int **F);
	~QA();
	int main();
	void initSpinConfig();
	void initBest();
	void SearchMaxDF();
	inline void selectSpin(int &d, int &a, int &b, int &p, int &q);
	inline double calcDeltaE(int d, int a, int b, int p, int q, double coefficient);
	inline void flipSpin(double coefficient);
	void printSpin();
	void printCost();
	inline void checkMin(int d);
private:
	random_device rnd;
	mt19937 mt;
	int ***spin;
	int bestCost=INT_MAX;
	int *bestPermutation;
	int *cost;
	int answer;
	double delta;
	int maxD = 0;
	int	maxF = 0;
	int **permutation;
	int **D, **F;
	int N;
	int trotterDim;
	double gamma;
	int annealingStep;
	int mcStep;
	double reducePara;
	double beta;//inverse T
};

