#include "common.h"
#include "SA.h"
#pragma once

class SearchParameter
{
public:
	SearchParameter(parameters p, int **D, int **F);
	~SearchParameter();
	void main(parameters &p, string modelName);
	int searchAnnealingStep(string modelName);
	int searchBeta(parameters p);
private:
	mt19937 mt;
	random_device rnd;
	int *x;
	double *result;
	int **D, **F;
	int maxD, maxF;
	const int trialsNum = 20;
	double alpha;
	int N;
	int R;
	double tInit;
	double tEnd;

	int getSuccess(int N, int *x, const double T);
	double meanT(double *result);
	vector<double> samplingDelta(int annealingStep);
	void searchMaxDF(int &maxD, int &maxF);
};