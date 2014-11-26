#pragma once
#include "../shadow_removal/decmps_proc.h"

struct DetectionCfg
{
	//enabled
	bool stepEnabled;

	//directory
	string srcDir;
	string paramSubdir;
	string resSubdir;

	//focused files
	vector<string> focusedPrefix;

	//step 1: mrf smooth
	bool enabled_mrf;
	int nLabels_mrf;
	int nCh_mrf;
	float colorsigma_mrfSTerm;
	float sWeight_mrfSTerm;

	//step 2: analysis
	DecmpsCfg dCfg;
	float guidInitMaxRatio_data, ratioSigma_data, whiteWeight_data, colorSigma_smooth, weight_smooth;

	//step 3: matting
	float shdwDegreeThres_matting;
	int maskErodeTimes_matting;
	int maskDilateTimes_matting;

	//debug
	bool debugImgsOutput;
	string debugImgsOutputDir;

	void Init();
	void InitFromXML(string cfgFile);
};

class DetectionProc
{
public:
	DetectionProc(void);
	~DetectionProc(void);

	void LoadCfg(string cfgFile);
	void DetectCastShadow();

	static DetectionCfg cfg;

	//abandon
	void SetFileDir(string fileDir);

private:

	cvi* DetectCastShadow(cvi* srcImg, cvi* param);

	cvi* MRFSmooth(cvi* srcImg, cvi* param, int nLabels = 64);
	cvi* DecmpsSmooth(cvi* srcImg, cvi* param, int nLabels = 64);

	cvi* MattingSmooth(cvi* srcImg, cvi* param, float shdwRatio = 0.5f, int erodeTimes = 5, int dilateTimes = 10);

private:

	//string m_fileDir;
	vector<string> m_imgNames;

	//string m_srcDir;
	//string m_resDir;
};

