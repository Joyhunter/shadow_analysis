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
	virtual ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20) = 0;
	virtual void ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info) = 0;

	ShdwImgInfo AnalysisShadowFrame(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN);
};

class RGBBiasMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20);
	void ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info);
};

class HLSGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20);
	void ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info);
};

class RGBBiasGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN = 20);
	void ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info);
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

