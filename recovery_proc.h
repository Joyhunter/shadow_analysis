#pragma once


class RecoverProc
{
public:
	RecoverProc(void);
	~RecoverProc(void);

	void SetFileDir(string fileDir);

	void Recover();

private:

	//smooth
	cvi* MRFSmooth(cvi* srcImg, cvi* param);
	cvi* DecmpsSmooth(cvi* srcImg, cvi* param);

	cvi* ImgRecoverNaive(cvi* srcImg, cvi* param);
	cvi* ImgRecoverPoisson(cvi* srcImg, cvi* param);

private:

	string m_fileDir;
	vector<string> m_imgNames;

	string m_ParamDir;
	string m_ResultDir;
	int m_patchRadius;
};

