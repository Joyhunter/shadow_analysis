#include "stdafx.h"
#include "shadow_analysis.h"
#include "gpm_analysis.h"

int _tmain(int argc, _TCHAR* argv[])
{
	
// 	ShdwAnlysisProc proc;
// 	proc.SetFileDir("..//dataset//newImg//indoor//");
// 	proc.ShdwAnlysis();

	GPMAnalysisProc gProc;
	gProc.SetFileDir("..//dataset//gpmAnalysis//1//");
	gProc.ShdwAnlysis();


	return 0;
}

