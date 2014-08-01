#pragma once

struct ShdwImgInfo
{
	int nPixels;
	//cvS avgBias, avgGain;
	vector<vector<int> > biasDstrbt, gainDstrbt;

	double asmptDist;
	
	ShdwImgInfo(int histN = 20):nPixels(0), biasDstrbt(3), gainDstrbt(3), asmptDist(0)
	{
		doF(k, 3)
		{
			biasDstrbt[k].resize(histN, 0);
			gainDstrbt[k].resize(histN, 0);
		}
	};
	void operator += (const ShdwImgInfo& info)
	{
		nPixels += info.nPixels;
		asmptDist += info.asmptDist;
		doFv(i, biasDstrbt) doFv(j, biasDstrbt[i])
		{
			biasDstrbt[i][j] += info.biasDstrbt[i][j];
		}
		doFv(i, gainDstrbt) doFv(j, gainDstrbt[i])
		{
			gainDstrbt[i][j] += info.gainDstrbt[i][j];
		}
	}
};

class Metric
{
public:
	virtual ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		int histN = 20, int patchRadius = 7) = 0;

	ShdwImgInfo AnalysisShadowBase(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN, 
		int patchRadius, vector<bool>& useGain, vector<bool>& useBias, int cSpaceType, int cSpaceType2);
	void ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info,
		vector<bool>& useGain, vector<bool>& useBias, int cSpaceType2);
};

class RGBGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20, int patchRadius = 7);
};

class HLSGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20, int patchRadius = 7);
};

class RGBBiasGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20, int patchRadius = 7);
};

class LABGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20, int patchRadius = 7);
};


class ShdwAnlysisProc
{
public:
	ShdwAnlysisProc(void);
	~ShdwAnlysisProc(void);

	void SetFileDir(string fileDir);

	void ShdwAnlysis();

private:

	cvi* VislzDstrbt(vector<vector<int> >& dstrbt, int nPixels);

private:

	string m_fileDir;
	vector<string> m_imgNames;

};

