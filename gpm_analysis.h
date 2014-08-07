#pragma once
#include "../shadow_removal/gpm_proc.h"

struct GPMParam
{
	GPMRange range;

	int downRatio, multiLevelRatio, extraLevelN; 
	int patchSize;

	int colorMode;
};

struct GPMStatis
{
	int matchN; //all matches
	int dssm; //different shadow, same material
	int sssm; //same shadow, same material
	int dsdm; //different shadow, different material
	int ssdm; //same shadow, different material
	GPMStatis():matchN(0), dssm(0), sssm(0), dsdm(0), ssdm(0){};
	void add(GPMStatis& v)
	{
		matchN += v.matchN; dssm += v.dssm; sssm += v.sssm; dsdm += v.dsdm; ssdm += v.ssdm;
	}
	void print()
	{
		cout<<"MatchN = "<<matchN<<".\n";
		cout<<"Dif shadow, same material\t= "<<dssm<<"("<<((matchN == 0)?0:(_f dssm*100/matchN))<<"%).\n";
		cout<<"Same shadow, same material\t= "<<sssm<<"("<<((matchN == 0)?0:(_f sssm*100/matchN))<<"%).\n";
		cout<<"Dif shadow, dif material\t= "<<dsdm<<"("<<((matchN == 0)?0:(_f dsdm*100/matchN))<<"%).\n";
		cout<<"Same shadow, dif material\t= "<<ssdm<<"("<<((matchN == 0)?0:(_f ssdm*100/matchN))<<"%).\n";
	}
};

class GPMAnalysisProc
{
public:
	GPMAnalysisProc(void);
	~GPMAnalysisProc(void);

	void SetFileDir(string fileDir);

	void ShdwAnlysis();

private:

	void RunGPMForAllImages();

	void GetStatistics();

	string GetCorrFileDir(string imgName);
	string GetShadowSegDir(string imgName);
	string GetMaterialSegDir(string imgName);

private:

	string m_fileDir;
	vector<string> m_imgNames;

	//param
	GPMParam param;
};

