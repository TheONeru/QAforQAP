#pragma once
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <time.h>
#include <omp.h>
#include <float.h>
#include <map>
#include<algorithm>
#include<numeric>
using namespace std;

typedef struct {
	int N;
	int trotterDim;
	double gamma;
	int annealingStep;
	int mcStep;
	double reducePara;
	double beta;
	int answer;
} parameters;

static void saveResult(string modelName, parameters p, vector<int> output, clock_t time) {

	int success = 0;
	long min_output = *min_element(output.begin(), output.end());
	long max_output = *max_element(output.begin(), output.end());
	double mean_output = 0;
	for (int i = 0; i < output.size(); i++) {
		if (output[i] <= p.answer) {
			success++;
		}
		mean_output += output[i];
	}
	mean_output /= output.size();

	cout << "success = " << success << endl;
	cout << "mean = " << mean_output << endl;
	cout << "best value = " << min_output << endl;
	cout << "worst value = " << max_output << endl;
	cout << "time = " << (double)time / CLOCKS_PER_SEC << endl;

	ofstream file("./result/" + modelName + ".txt");
	file << "success = " << success << endl;
	file << "mean = " << mean_output << endl;
	file << "best value = " << min_output << endl;
	file << "worst value = " << max_output << endl;
	file << "time = " << (double)time / CLOCKS_PER_SEC << endl;
	file << "annelingStep = " << p.annealingStep << endl;
	file << "beta = " << p.beta << endl;
	file << "trotterDim = " << p.trotterDim << endl;
	file.close();
}

static parameters setParameter(string modelName) {
	parameters p;
	if (modelName == "rou15") {
		p.answer = 354210;
		p.N = 15;
	}
	else if (modelName == "had12") {
		p.answer = 1652;
		p.N = 12;
	}
	else if (modelName == "chr20c") {
		p.answer = 14142;
		p.N = 20;
	}
	else if (modelName == "had20") {
		p.answer = 6922;
		p.N = 20;
	}
	else if (modelName == "tai25a") {
		p.answer = 1167256;
		p.N = 25;
	}
	else if (modelName == "tai25b") {
		p.answer = 344355646;
		p.N = 25;
	}
	else if (modelName == "bur26a") {
		p.answer = 5426670;
		p.N = 26;
	}
	else if (modelName == "kra30b") {
		p.answer = 91420;
		p.N = 30;
	}
	else if (modelName == "tai30a") {
		p.answer = 1818146;
		p.N = 30;
	}
	else if (modelName == "nug30") {
		p.answer = 6124;
		p.N = 30;
	}
	else if (modelName == "tai30b") {
		p.answer = 637117113;
		p.N = 30;
	}
	else if (modelName == "tai40b") {
		p.answer = 637250948;
		p.N = 40;
	}
	else if (modelName == "tai50a") {
		p.answer = 4938796;
		p.N = 50;
	}
	else if (modelName == "tai50b") {
		p.answer = 458821517;
		p.N = 50;
	}
	else if (modelName == "tai80a") {
		p.answer = 13499184;
		p.N = 80;
	}
	p.trotterDim = p.N;//調べられていないが，十分大きければOK，今回はSAに合わせるために
	p.mcStep = 10*p.N*p.N;//SAにおけるR(平衡状態作成のループ数)と同じに
	p.gamma = 1;//調べられていないが，とりあえず1.0
	p.reducePara = 0.95;//SAの減衰係数を参考に
	/*文献を調べられていない，参考書には十分に大きくとだけ(betaは不変)，
	  最終アニーリングステップで,仮にある状態からの差のサンプリングが正規分布だと仮定した時の99.73%が横磁場の影響を受けるような設定
	*/
	p.beta = 10000;
	p.annealingStep = 2000;//SAと計算量を同じに，後で計算もしくはファイルで検索
	return p;
}

static void getMapdatas(int N, string modelName, int**D, int **F) {
	FILE *file;
	string path = "constants/" + modelName + ".txt";
	fopen_s(&file, path.c_str(), "r");
	if (!file) {
		cout << "Unkwon Model name" + modelName << endl;
		exit(1);
	}
	string str;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fscanf_s(file, "%d", &D[i][j]);
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fscanf_s(file, "%d", &F[i][j]);
		}
	}
	fclose(file);
}