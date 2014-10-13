#pragma once

class DetectionProc
{
public:
	DetectionProc(void);
	~DetectionProc(void);

	void SetFileDir(string fileDir);

	void DetectCastShadow();

private:

	cvi* DetectCastShadow(cvi* srcImg, cvi* param);

	cvi* MRFSmooth(cvi* srcImg, cvi* param, int nLabels = 64);
	cvi* DecmpsSmooth(cvi* srcImg, cvi* param, int nLabels = 64);

private:

	string m_fileDir;
	vector<string> m_imgNames;

	string m_srcDir;
	string m_resDir;
};

