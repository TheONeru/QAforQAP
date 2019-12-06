#pragma once
#include "common.h"
class QA
{
public:
	QA(parameters p, int **D, int **F);
	~QA();
	int main();
	void initPermutation();
	void getSpinConfig();
	void initBest();
	void SearchMaxDF();
	inline void selectSpin(int &d, int &a, int &b, int &p, int &q);
	inline void selectPositon(int &a, int &b);
	inline double calcDeltaE(int d, int a, int b, int p, int q, double coefficient);
	inline void flipSpin(double coefficient);
	inline void sa(double invT);
	inline void checkMin(int d);
	int SAQA();
	int expMain(string modelName);

private:
	random_device rnd;
	mt19937 mt;
	int ***spin;
	int bestCost=INT_MAX;
	double initT;
	double reduceT;
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
	double reduceGamma;
	double beta;//inverse T
	int acceptMinus = 0;
	int acceptMetro = 0;
	double tmpDoubleE = 0;
	double tmpSougo = 0;
};

