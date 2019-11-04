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
	p.trotterDim = p.N;//���ׂ��Ă��Ȃ����C�\���傫�����OK�C�����SA�ɍ��킹�邽�߂�
	p.mcStep = 10*p.N*p.N;//SA�ɂ�����R(���t��ԍ쐬�̃��[�v��)�Ɠ�����
	p.gamma = 1;//���ׂ��Ă��Ȃ����C�Ƃ肠����1.0
	p.reducePara = 0.95;//SA�̌����W�����Q�l��
	/*�����𒲂ׂ��Ă��Ȃ��C�Q�l���ɂ͏\���ɑ傫���Ƃ���(beta�͕s��)�C
	  �ŏI�A�j�[�����O�X�e�b�v��,���ɂ����Ԃ���̍��̃T���v�����O�����K���z���Ɖ��肵������99.73%��������̉e�����󂯂�悤�Ȑݒ�
	*/
	p.beta = 10000;
	p.annealingStep = 2000;//SA�ƌv�Z�ʂ𓯂��ɁC��Ōv�Z�������̓t�@�C���Ō���
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