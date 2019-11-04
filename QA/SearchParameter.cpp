#include "SearchParameter.h"

SearchParameter::SearchParameter(parameters p, int **D, int **F)
{
	mt.seed(rnd());
	N = p.N;
	R = p.mcStep;
	alpha = p.reducePara;
	this->D = new int*[N];
	this->F = new int*[N];
	for (int i = 0; i < N; i++) {
		this->D[i] = new int[N];
		this->F[i] = new int[N];
		memcpy(this->D[i], D[i], sizeof(int)*N);
		memcpy(this->F[i], F[i], sizeof(int)*N);
	}
	result = new double[trialsNum];
	x = new int[N];
}

SearchParameter::~SearchParameter()
{
	for (int i = 0, N = this->N; i < N; i++) {
		delete[] D[i], F[i];
	}
	delete[] F, D, x, result;
}

void SearchParameter::main(parameters &p, string modelName) {
	p.annealingStep = searchAnnealingStep(modelName);
	cout << "初期温度 = " << tInit << endl;
	cout << "終了温度 = " << tEnd << endl;
	cout << "アニーリングステップ = " << p.annealingStep << endl;
	p.beta = searchBeta(p);
	cout << "目標値 平均+3σ(正規分布99.73%)での逆温度 = " << p.beta << endl;
}

int SearchParameter::searchAnnealingStep(string modelName) {
	ifstream parameterFile("./constants/parameters/" + modelName + ".csv");
	if (parameterFile.is_open()) {
		parameterFile>>tInit;
		parameterFile>>tEnd;
		parameterFile.close();
		return (int)N*(log(tEnd / tInit) / log(alpha));
	}
	parameterFile.close();
	ofstream outFile("./constants/parameters/" + modelName + ".csv");

	int N = this->N;
	double alpha = this->alpha;
	int p_init, p_end;
	int success = 0;//状態が改善した回数
	
    //受け入れ確立
	p_init = (int)R * 0.95;
	p_end = (int)R * 0.005;

	int *x = this->x;
	double *result = this->result;
	//p_initの温度探索
	for (int i = 0, n = trialsNum; i < n; i++) {
		createInitState(N, x, mt);
		double T = LONG_MAX;
		T *= 2;
		for (long changeT = LONG_MAX; changeT > 2; changeT >>= 1) {
			success = getSuccess(N, x, T);
			if (success <= p_init) {
				T += changeT;
			}
			else {
				T -= changeT;
			}
		}
		result[i] = T;
	}
	tInit = meanT(result);

	double tLim = 0.1;
	bool flg = false;
	createInitState(N, x, mt);
	//p_endの温度探索
	for (int i = 0, n = trialsNum; i < n; i++) {
		double T;
		for (T = tInit; T > tLim; T *= alpha) {
			success = getSuccess(N, x, T);
			if (success <= p_end) {
				break;
			}
		}
		result[i] = T;

	}
	tEnd = meanT(result);
	outFile << tInit << " " << tEnd << endl;
	return (int)N*(log(tEnd / tInit) / log(alpha));
}

int SearchParameter::searchBeta(parameters p) {

	vector<double> sample = samplingDelta(p.annealingStep);

	//検証用
	/*ofstream sampleFile("./exp/samplingData.csv");
	for (int i = 0; i < sample.size(); i++) {
		sampleFile << sample[i] << endl;
	}
	sampleFile.close();*/

	double sum = 0;
	double sum2 = 0;//標準偏差用
	for (int i = 0, n=sample.size(); i < n; i++) {
		sum += sample[i];
		sum2 += sample[i] * sample[i];
	}
	double mean = sum / sample.size();
	double dev = sqrt(sum2 / sample.size());
	cout << "mean=" << mean <<endl;
	double magnification = -1*N*((double)2/250) + (double)29/10;
	double target = mean + 2.4 * dev;
	cout << "目標値 = " << target << endl;
	double betaAns = LONG_MAX;
	double tmpA = pow(p.reducePara, p.annealingStep) / p.trotterDim;
	for (double b = LONG_MAX; b > 1; b /= 2) {
		double c = (1 / betaAns)*log(tanh(betaAns*tmpA));
		if (-4 * c < target) {
			betaAns -= b;
		}
		else {
			betaAns += b;
		}
	}
	return betaAns;
}

inline int SearchParameter::getSuccess(int N, int *x, double T) {
	int select1, select2;
	int R = this->R;
	double T_inv = 1 / T;
	int **D = this->D;
	int **F = this->F;
	int cnt = 0;
	for (int i = 0; i < R; i++) {
		createNewState(N, select1, select2, mt);
		int d = calcDelta(N, x, D, F, select1, select2);
		int buff;
		uniform_real_distribution<double>rand(0, 1);
		if (d <= 0) {
			cnt++;
			buff = x[select1];
			x[select1] = x[select2];
			x[select2] = buff;
		}
		else if (rand(mt) <= exp(-d * T_inv)) {
			cnt++;
			buff = x[select1];
			x[select1] = x[select2];
			x[select2] = buff;
		}
	}
	return cnt;
}

double SearchParameter::meanT(double *result) {
	double ave = 0;
	double min = ~(1 << 31), max = 0;
	for (int i = 0, n = trialsNum; i < n; i++) {
		ave += result[i];
		if (min > result[i]) {
			min = result[i];
		}
		if (max < result[i]) {
			max = result[i];
		}
	}
	ave -= (min + max);
	ave /= trialsNum;
	return ave;
}

vector<double> SearchParameter::samplingDelta(int annealingStep) {
	vector<double> samplingData;
	int select1, select2;
	int R = this->R;
	int **D = this->D;
	int **F = this->F;
	long f = 0;
	double tEnd = this->tEnd;
	double alpha = this->alpha;
	double T = tInit;
	int maxD=this->maxD;
	int maxF=this->maxF;
	createInitState(N, x, mt);
	initCost(N, x, D, F, f);
	searchMaxDF(maxD, maxF);
	for (int i = 0; i < N*1000; i++) {
		createInitState(N, x, mt);
		createNewState(N, select1, select2, mt);
		int d = calcDelta(N, x, D, F, select1, select2);
		samplingData.push_back((double)d / (double)(maxD*maxF));
	}
	return samplingData;
}

void SearchParameter::searchMaxDF(int &maxD, int &maxF) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (maxD < D[i][j]) {
				maxD = D[i][j];
			}
			if (maxF < F[i][j]) {
				maxF = F[i][j];
			}
		}
	}
}