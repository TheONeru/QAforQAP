#include "QA.h"

QA::QA(parameters p, int **D, int **F)
{
	mt.seed(rnd());
	N = p.N;
	trotterDim = p.trotterDim;
	gamma = p.gamma;
	annealingStep = p.annealingStep;
	mcStep = p.mcStep;
	reducePara = p.reducePara;
	beta = p.beta;
	answer = p.answer;
	//spin array allocate memory
	spin = new int**[N];
	permutation = new int*[trotterDim];
	for (int t = 0, trotterDim = this->trotterDim; t < trotterDim; t++) {
		spin[t] = new int*[N];
		for (int l = 0; l < N; l++) {
			spin[t][l] = new int[N];
		}
		permutation[t] = new int[N];
	}
	this->D = new int*[N];
	this->F = new int*[N];
	for (int i = 0; i < N; i++) {
		this->D[i] = new int[N];
		this->F[i] = new int[N];
		memcpy(this->D[i], D[i], sizeof(int)*p.N);
		memcpy(this->F[i], F[i], sizeof(int)*p.N);

	}
	this->cost = new int[trotterDim];
	this->bestPermutation = new int[N];
	
}

QA::~QA()
{
	//free memory
	for (int t = 0, trotterDim = this->trotterDim; t < trotterDim; t++) {
		for (int l = 0, N = this->N; l < N; l++) {
			delete[] spin[t][l];
		}
		delete spin[t];
		delete permutation[t];
	}
	for (int t = 0, N = this->N; t < N; t++) {
		delete D[t],F[t];
	}
	delete[] spin, cost, F, D, permutation;
}

int QA::main() {	
	double gamma = this->gamma;
	double reducePara = this->reducePara;
	int answer = this->answer;
	initSpin();
	initBest();
	SearchMaxDF();
	for (int s = 0, annealingStep = this->annealingStep; s < annealingStep; s++) {
		double coefficient = (1 / beta)*log(tanh(beta*gamma / trotterDim));
		for (int m = 0, mcStep = this->mcStep; m < mcStep; m++) {
			flipSpin(coefficient);
		}
		gamma *= reducePara;
		if (bestCost <= answer) {
			break;
		}
	}
	for (int i = 0, N=this->N; i < N; i++) {
		cout << bestPermutation[i] << " ";
	}
	cout << endl;
	cout << bestCost << endl;
	return bestCost;
}

//初期スピン配置の決定, 配列の次元は[Trotter, location(vertical), factory(horizon)]
/*-
下の例のように一つだけ1に残りは-1にする(1,-1にすることで反転を容易にしている)
-1,-1,1
1,-1,-1
-1,1,-1
*/
void QA::initSpin() {
	int ***spin = this->spin;

	//spin変更候補作成
	vector<int> spinCandidate(N);
	iota(spinCandidate.begin(), spinCandidate.end(), 0);

	//初期スピン配置作成
	for (int t = 0, trotterDim = this->trotterDim; t < trotterDim; t++) {
		shuffle(spinCandidate.begin(), spinCandidate.end(), mt);
		for (int l = 0, N = this->N; l < N; l++) {
			for (int f = 0; f < N; f++) {
				spin[t][l][f] = -1;
			}
			int s_ = spinCandidate[l];
			spin[t][l][s_] = 1;
		}
	}

	//データ構造変更，計算上順列の方が計算しやすいため順列に変更
	for (int t = 0, trotterDim = this->trotterDim; t < trotterDim; t++) {
		for (int l = 0, N = this->N; l < N; l++) {
			for (int f = 0; f < N; f++) {
				if (spin[t][l][f] == 1) {
					permutation[t][l] = f;
					break;
				}
			}
		}
	}	

}

void QA::initBest() {
	int **permutation = this->permutation;
	for (int t = 0, trotterDim = this->trotterDim; t < trotterDim; t++) {
		int cost = 0;
		for (int l1 = 0, N = this->N; l1 < N; l1++) {
			for (int l2 = 0; l2 < N; l2++) {
				cost += (long)D[l1][l2] * F[permutation[t][l1]][permutation[t][l2]];
			}
		}
		this->cost[t] = cost;
		if (bestCost > cost) {
			bestCost = cost;
			memcpy(bestPermutation, permutation[t], sizeof(int)*N);
		}
	}
}

//数式には書かれていないが，QAPでは差分が大きくなりすぎて規格化しないと横磁場が掛かりにくく収束しにくいため最大のD，F要素をみつける
void QA::SearchMaxDF() {
	int &maxD = this->maxD;
	int &maxF = this->maxF;
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

//ひっくり返すスピンの選択
/*
d次元目のa,bがそれぞれ行の0,1に対応した時にa行目の1がある列をp, b行目の1がある列をq
   p,q
  -----
a| 1,-1
b| -1,1
*/
inline void QA::selectSpin(int &d, int &a, int &b, int &p, int &q) {
	int ***spin = this->spin;
	uniform_int_distribution<> random_dim(0, trotterDim - 1);
	d = random_dim(mt);
	uniform_int_distribution<> random_location(0, N-1);
	a = random_location(mt);
	b = random_location(mt);
	while (a == b) {
		b = random_location(mt);
	}
	//二分木探索にしたらもっと早くなりそうだけど，問題サイズが比較的小さいからあまり変わらなさそうでもある
	for (int f = 0, N = this->N; f < N; f++) {
		if (spin[d][a][f] == 1) {
			p = f;
			break;
		}
	}
	for (int f = 0, N = this->N; f < N; f++) {
		if (spin[d][b][f] == 1) {
			q = f;
			break;
		}
	}
}

//目的関数全体の差分
inline double QA::calcDeltaE(int d, int a, int b, int p, int q, double coefficient) {
	int **permutation = this->permutation;
	int **D = this->D;
	int **F = this->F;
	delta = 0;
	for (int i = 0, N = this->N; i < N; i++) {
		delta += (D[i][a] - D[i][b])*(F[q][permutation[d][i]] - F[p][permutation[d][i]]);
		delta += (D[a][i] - D[b][i])*(F[permutation[d][i]][q] - F[permutation[d][i]][p]);
	}
	delta += (D[a][b] + D[b][a] - D[a][a] - D[b][b])*(F[p][q] + F[q][p] - F[p][p] - F[q][q]);
	delta /= (maxD*maxF);

	//負の余剰は問題があるので正にしてから余剰にしている
	//それぞれの一変数の変化による差分, 第2項
	double delta1 = spin[d][a][p] * (spin[(d + trotterDim - 1) % trotterDim][a][p] + spin[(d + 1) % trotterDim][a][p]);
	double delta2 = spin[d][a][q] * (spin[(d + trotterDim - 1) % trotterDim][a][q] + spin[(d + 1) % trotterDim][a][q]);
	double delta3 = spin[d][b][p] * (spin[(d + trotterDim - 1) % trotterDim][b][p] + spin[(d + 1) % trotterDim][b][p]);
	double delta4 = spin[d][b][q] * (spin[(d + trotterDim - 1) % trotterDim][b][q] + spin[(d + 1) % trotterDim][b][q]);
	
	double deltaE=0;
	//tmpSougo = (delta1 + delta2 + delta3 + delta4);
	deltaE = (double)delta / trotterDim + coefficient*(delta1 + delta2 + delta3 + delta4);
	return deltaE;
}

inline void QA::flipSpin(double coefficient) {
	int ***spin = this->spin;
	int d, a, b, p, q;
	selectSpin(d, a, b, p, q);
	double delta_E = calcDeltaE(d, a, b, p, q, coefficient);
	tmpDoubleE = delta_E;
	uniform_real_distribution<> r(0, 1);
	if (delta_E <= 0) {
		spin[d][a][p] *= -1;
		spin[d][a][q] *= -1;
		spin[d][b][p] *= -1;
		spin[d][b][q] *= -1;
		permutation[d][a] = q;
		permutation[d][b] = p;
		cost[d] += delta* (maxD*maxF);
	}
	else if (r(mt)<exp(-beta*delta_E)) {
		spin[d][a][p] *= -1;
		spin[d][a][q] *= -1;
		spin[d][b][p] *= -1;
		spin[d][b][q] *= -1;
		permutation[d][a] = q;
		permutation[d][b] = p;
		cost[d] += delta * (maxD*maxF);
	}
	checkMin(d);
}

inline void QA::checkMin(int d) {
	if (bestCost > cost[d]) {
		bestCost = cost[d];
		memcpy(bestPermutation, permutation[d], sizeof(int)*N);
	}
}

int QA::expMain(string modelName) {
	ofstream outputFile("./exp/" + modelName + ".csv");
	ofstream acceptFile("./exp/" + modelName + "Accept.csv");
	//ofstream deltaFile("./exp/" + modelName + "Delta.csv");
	//ofstream sougoFile("./exp/" + modelName + "Sougo.csv");
	double gamma = this->gamma;
	double reducePara = this->reducePara;
	int answer = this->answer;
	initSpin();
	initBest();
	SearchMaxDF();
	for (int s = 0, annealingStep = this->annealingStep; s < annealingStep; s++) {
		double coefficient = (1 / beta)*log(tanh(beta*gamma / trotterDim));
		cout << s << "annealler step" << endl;
		for (int m = 0, mcStep = this->mcStep; m < mcStep; m++) {
			flipSpin(coefficient);
			if (m % 2500==0) {
				outputFile << s * mcStep + m << "," << cost[0] << endl;
				acceptFile << s * mcStep + m << "," << tmpDoubleE <<endl;
				//deltaFile << s * mcStep + m << "," << delta << endl;
				//sougoFile << s * mcStep + m << "," << tmpSougo << endl;
			}
		}
		gamma *= reducePara;
	}
	for (int i = 0, N = this->N; i < N; i++) {
		cout << bestPermutation[i] << " ";
	}
	cout << endl;
	cout << bestCost << endl;
	return bestCost;
}
