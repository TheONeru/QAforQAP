#include "QAManager.h"

QAManager::QAManager(string modelName, double tInit, double tEnd)
{
	this->modelName = modelName;
	p = setParameter(modelName, tInit, tEnd);
	this->N = p.N;
	int N = p.N;
	D = new int*[N];
	F = new int*[N];
	for (int i = 0; i < N; i++) {
		D[i] = new int[N];
		F[i] = new int[N];
	}
	getMapdatas(N, modelName, D, F);
}

QAManager::~QAManager()
{
	for (int t = 0, N = this->N; t < N; t++) {
		delete D[t], F[t];
	}
	delete[] F, D;
}

void QAManager::exe() {
	int trialNum = 100;
	vector<int> solution(trialNum, 0);
	clock_t start = clock();
#pragma omp parallel for
	for (int i = 0; i < trialNum; i++) {
		QA q(p, D, F);
		solution[i] = q.main();
	}
	clock_t end = clock();
	saveResult(modelName, p, solution, (end - start));
}
