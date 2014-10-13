#pragma once

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

	void SetFileDir(string fileDir);

	//void Recover();
	//void Recover2();
	void Recover3();

private:

	//naive direct compute result
	cvi* ImgRecoverNaive(cvi* srcImg, cvi* param, int patchSize = 7);

	//synthesis
	void GenerateMaskFromParam(IN cvi* param, OUT cvi* holeMask, OUT cvi* legalMask, float thres = 0.8f, float bound = 0.05f);



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

	void GetScript();

private:

	string m_fileDir;
	vector<string> m_imgNames;

	string m_ParamDir;
	string m_ResultDir;
	int m_patchRadius;
};

