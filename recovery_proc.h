#pragma once
#include "synthesis_proc.h"

struct RecoveryCfg
{
	//enabled
	bool stepEnabled;

	//directory
	string srcDir;
	string paramSubdir;
	string resSubdir;

	//focused files
	vector<string> focusedPrefix;

	//step 1: compute naive res
	string naiveResImgPostfix;
	int patchsize_naiveRecovery;

	//step 2: synthesis
	//mask generation
	float shdwDegreeThres_maskGenerate;
	float shdwBoundDilateTimes_maskGenerate;
	string holeMaskImgPostfix, legalMaskImgPostfix;
	//synthesis
	//SynthesisCfg synCfg;
	bool gtAsGuidance;
	string synResImgPostfix;

	//step 3: local color correction
	int patchRadius_localColorCorrection;
	float correctionStepRatio_localColorCorrection;
	string corResImgPostfix;
	//pyramid
	int pyramidExtraLevels_pyramid;
	float pyramidResizeRatio_pyramid;

	//origial
	int patchRadius_old;

	void Init();
	void InitFromXML(string cfgFile);
};

struct ImgPrmd
{
	vector<cvi*> imgs;

	ImgPrmd(cvi* img, int levels, float ratio);
	cvi* Flatten();
	void Release();
};

class RecoverProc
{
public:
	RecoverProc(void);
	~RecoverProc(void);

	void LoadCfg(string cfgFile);

	//void Recover();
	//void Recover2();
	void Recover3();

	//naive direct compute result
	static cvi* ImgRecoverNaive(cvi* srcImg, cvi* param, int patchSize = 7);

	static RecoveryCfg cfg;

private:


	//synthesis
	void GenerateMaskFromParam(IN cvi* param, OUT cvi* holeMask, OUT cvi* legalMask, float thres = 0.8f, float bound = 0.05f);

	cvi* LocalColorCorrection(cvi* naiveRes, cvi* synRes, cvi* holeMask);
	cvi* LocalColorCorrectionSingleLevel(cvi* naiveRes, cvi* synRes, cvi* holeMask, int patchRadius = 10);

	cvi* GetCorCpltPixels(cvi* naiveRes, cvi* synRes);



	//smooth
	//cvi* MRFSmooth(cvi* srcImg, cvi* param, int nLabels = 64);
	//cvi* DecmpsSmooth(cvi* srcImg, cvi* param, int nLabels = 64);
	cvi* BoundarySmooth(cvi* srcImg, cvi* smoothParam, cvi* &shadowMask, float Lthres = 200.f, float r1 = 5.f, float r2 = 10.f);

	cvi* VislzBound(cvi* src, cvi* boundMask);

	cvi* ImgRecoverPoisson(cvi* srcImg, cvi* param);

	void PoissonSmooth(cvi* img, cvi* mask);

	//whole image, optimize best param
	cvi* GetBestParam(IN cvi* src, IN cvi* mask, IN cvi* cpImg, OUT cvS& param);
	cvi* GetBestParamLapPymd(IN cvi* src, IN cvi* mask, IN cvi* cpImg, OUT cvS& param);

	//abandon
	void SetFileDir(string fileDir);
	void GetScript();

private:


	vector<string> m_imgNames;

};

