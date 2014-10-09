#include "stdafx.h"
#include "shadow_analysis.h"
#include "gpm_analysis.h"
#include "recovery_proc.h"
#include "synthesis_proc.h"

int _tmain(int argc, _TCHAR* argv[])
{
	
// 	cvi* test = cvlic("316__v2.png");
// 	
// 	cvi* src2 = cvci81(test);
// 
// 	doFcvi(test, i, j)
// 	{
// 		cvs20(src2, i, j, cvg20(test, i, j));
// 	}
// 
// 	cvsi("316__v2.png", src2);
// 	return 0;


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

	SynthesisProc sProc;
	sProc.Test();

// 	RecoverProc rProc;
// 	rProc.SetFileDir("..//dataset//gpmAnalysis//3//");
// 	//rProc.Recover();
// 	rProc.Recover2();


	return 0;
}


// 	proc.AddFileDir("..//dataset//myData//indoor//1//");
// 	proc.AddFileDir("..//dataset//otherData//indoor//");
// 	proc.AddFileDir("..//dataset//otherData//outdoor//");
// 	proc.AddFileDir("..//dataset//myData//indoor//2//");
// 	proc.AddFileDir("..//dataset//myData//outdoor//1//");
// 	proc.AddFileDir("..//dataset//myData//outdoor//2//");