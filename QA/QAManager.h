#pragma once
#include "QA.h"
#include "common.h"

class QAManager
{
public:
	QAManager(string modelName, double tInit, double tEnd);
	~QAManager();
	void exe();
private:
	string modelName;
	parameters p;
	int N;
	int **D;
	int **F;
};

