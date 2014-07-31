#include "stdafx.h"
#include "shadow_analysis.h"

int _tmain(int argc, _TCHAR* argv[])
{
	
	ShdwAnlysisProc proc;
	proc.SetFileDir("..//dataset//newImg//outdoor//");
	proc.ShdwAnlysis();

	return 0;
}

