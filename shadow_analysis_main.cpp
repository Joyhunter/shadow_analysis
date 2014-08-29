#include "stdafx.h"
#include "shadow_analysis.h"
#include "gpm_analysis.h"
#include "recovery_proc.h"

int _tmain(int argc, _TCHAR* argv[])
{
	
//     ShdwAnlysisProc proc;
//  	proc.AddFileDir(".//");
//  	proc.BunchShdwAnalysis();
// 	//proc.GetStatistics();
// 	return 0;

// 	GPMAnalysisProc gProc;
// 	//gProc.SetFileDir("..//dataset//gpmAnalysis//2//");
// 	gProc.SetFileDir(".//");
// 	gProc.ShdwAnlysis();
// 	//gProc.VoteAnalysis();
// 	return 0;

	RecoverProc rProc;
	rProc.SetFileDir("..//dataset//gpmAnalysis//3//");
	rProc.Recover();


	return 0;
}


// 	proc.AddFileDir("..//dataset//myData//indoor//1//");
// 	proc.AddFileDir("..//dataset//otherData//indoor//");
// 	proc.AddFileDir("..//dataset//otherData//outdoor//");
// 	proc.AddFileDir("..//dataset//myData//indoor//2//");
// 	proc.AddFileDir("..//dataset//myData//outdoor//1//");
// 	proc.AddFileDir("..//dataset//myData//outdoor//2//");