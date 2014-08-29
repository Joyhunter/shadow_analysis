#pragma once

struct ShdwImgInfo
{
	int nPixels, nImages;
	//cvS avgBias, avgGain;
	vector<vector<int> > biasDstrbt, gainDstrbt;

	vector<bool> useBias, useGain;
	vector<double> biasMin, biasMax, gainMin, gainMax;

	double asmptDist;
	
	ShdwImgInfo(int histN = 20):nPixels(0), nImages(0), biasDstrbt(3), gainDstrbt(3), asmptDist(0)
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
		useBias = info.useBias;
		useGain = info.useGain;
		nImages += info.nImages;
	}

	int computeBiasIdx(double bias, int histN, bool useBias, bool useGain);
	double computeBiasV(int biasIdx, int histN, bool useBias, bool useGain);
	int computeGainIdx(double gain, int histN, bool useBias, bool useGain);
	double computeGainV(int gainIdx, int histN, bool useBias, bool useGain);

	void computeMinMax(float minMaxRatio = 0.1f);
	void outputMinMax();
};


class ImgConverter
{
public:
	virtual void convert(const cvi* src, cvi* dst) = 0;
	virtual void invConvert(const cvi* src, cvi* dst) = 0;
};

class RGBRatioConverter : public ImgConverter
{
public:
	void convert(const cvi* src, cvi* dst);
	void invConvert(const cvi* src, cvi* dst);
};


class Metric
{
public:
	virtual ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		OUT vector<double>& gbVs, int histN = 20, int patchRadius = 7) = 0;

	ShdwImgInfo AnalysisShadowBase(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN, 
		int patchRadius, vector<bool>& useGain, vector<bool>& useBias, int cSpaceType, int cSpaceType2, 
		OUT vector<double>& gbVs, ImgConverter* converter = NULL);
	void ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info,
		vector<bool>& useGain, vector<bool>& useBias, int cSpaceType2, OUT vector<int>& idxs, OUT vector<double>& gbV,
		ImgConverter* converter = NULL);
};

class RGBGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		OUT vector<double>& gbVs, int histN = 20, int patchRadius = 7);
};

class HLSGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		OUT vector<double>& gbVs, int histN = 20, int patchRadius = 7);
};

class RGBBiasGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		OUT vector<double>& gbVs, int histN = 20, int patchRadius = 7);
};

class LABGainMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		OUT vector<double>& gbVs, int histN = 20, int patchRadius = 7);
};

//kalyan's advise 2
class RGBRatioMetric : public Metric
{
public:
	ShdwImgInfo AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
		OUT vector<double>& gbVs, int histN = 20, int patchRadius = 7);
};

class ShdwAnlysisProc
{
public:
	ShdwAnlysisProc(void);
	~ShdwAnlysisProc(void);

	void AddFileDir(string fileDir);

	void BunchShdwAnalysis();
	void GetStatistics();

	static cvi* VislzVector(int width, int height, vector<double>& vs, int offset, double minV, double maxV);
	static void SaveVector(int width, int height, vector<double>& vs, int offset, string fileDir);

	static cvi* VislzFCvi(int width, int height, cvi* src, int ch, double minV, double maxV);
	static void SaveFCvi(int width, int height, cvi* src, int ch, string fileDir);


private:

	ShdwImgInfo ShdwAnlysis(string fileDir, IN int histN, IN float resizeRatio, IN int patchRadius,
		IN float minMaxRatio, IN Metric* metric);

	cvi* VislzDstrbt(vector<vector<int> >& dstrbt, int nPixels);


private:

	vector<string> m_fileDirs;
	vector<string> m_imgNames;

};

