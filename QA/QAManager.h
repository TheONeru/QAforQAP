#pragma once
#include "QA.h"
#include "common.h"
#include "SearchParameter.h"

class QAManager
{
public:
	QAManager(string modelName);
	~QAManager();
	void exe();
	void expExe();
private:
	string modelName;
	parameters p;
	int N;
	int **D;
	int **F;
};

