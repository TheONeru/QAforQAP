#include "common.h"
#pragma once

static void initCost(int N, int* x, int** const D, int** const F, long &f) {
	f = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f += (long)D[i][j] * F[x[i]][x[j]];
		}
	}
}

static inline void createInitState(int N, int *x, mt19937 &mt) {
	vector<int> randomPermutation(N);
	for (int j = 0; j < N; j++) {
		randomPermutation[j] = j;
	}
	for (int j = 0; j < N; j++) {
		uniform_int_distribution<> initRandom(0, randomPermutation.size() - 1);
		int r = initRandom(mt);
		x[j] = randomPermutation[r];
		randomPermutation.erase(randomPermutation.begin() + r);
	}
}

static inline void createNewState(int N, int &select1, int &select2, mt19937 &mt) {
	uniform_int_distribution<> rand(0, N - 1);
	select1 = rand(mt);
	select2 = rand(mt);
	while (select1 == select2) {
		select1 = rand(mt);
		select2 = rand(mt);
	}
}

static inline int calcDelta(int N, int* const x, int** const D, int** const F, int P1, int P2) {
	int delta = 0;
	for (int i = 0; i < N; i++) {
		delta += (D[i][P1] - D[i][P2])*(F[x[P2]][x[i]] - F[x[P1]][x[i]]);
		delta += (D[P1][i] - D[P2][i])*(F[x[i]][x[P2]] - F[x[i]][x[P1]]);
	}
	delta += (D[P2][P1] + D[P1][P2] - D[P2][P2] - D[P1][P1])
		*(F[x[P1]][x[P2]] + F[x[P2]][x[P1]] - F[x[P1]][x[P1]] - F[x[P2]][x[P2]]);
	return delta;
}
