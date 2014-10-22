#pragma once

struct InputImageData
{
	cvi* src;
	cvi* holeMask;
	cvi* legalMask;
	cvi* guideImg;
	InputImageData(cvi* src, cvi* holeMask, cvi* legalMask, cvi* guideImg):src(src), holeMask(holeMask),
		legalMask(legalMask), guideImg(guideImg){};
};

class DenseCorrSyn;

struct SynthesisCfg
{
	int patchSize; // 7

	//synthesis
	float pymResizeRatio; // 0.5
	int pymLevels; // 5
	int cpltItrN; // 2
	float poissonAlpha; // 0.1

	//gpm
	int gpmItrN; // 2
	bool hFlipEnabled; // true
	bool vFlipEnabled; // true
	bool scaleRotateEnabled; // false
	bool gainBiasEnabled; // false
	bool randomSearchEnabled; // true
	float randomSearchScaleFactor; // 0.5
	float randomSearchRadius; // 0.1 0.2

	//feature
	float guideImgWeight; // 0.0

	void Init();
};

class SynthesisProc
{

public:

	SynthesisProc(void);
	~SynthesisProc(void);

	void Synthesis(cvi* srcImg, cvi* holeMask, cvi* resImg, cvi* legalMask = NULL, cvi* guideImg = NULL);

	void Test();

	static SynthesisCfg cfg;

private:

	void VoteToCompletion(DenseCorrSyn* dsCor, InputImageData& imageData, cvi* resImg, int patchSize, float gradientAlpha = 1);	

	void PossionSolve(cvi* srcImg, cvi* holeImg, cvi* resImg, cvi* resConstrain, 
		cvi* gradientConstrainV, cvi* gradientConstrainH, float gradientAlpha = 1);


};


//--------------------- Interval & GPMRange ------------------------------------

class IntervalSyn
{
public:
	float min, max;
	IntervalSyn(float min = 0, float max = 0);
	float RandValue();
	friend ostream& operator << (ostream& out, IntervalSyn& val){
		out<<val.min<<" "<<val.max;
		return out;
	}
};

struct GPMSynRange
{
	IntervalSyn m_scaleItrl, m_rotateItrl;
	vector<IntervalSyn> m_biasItrl, m_gainItrl;
	bool hrEnable, vrEnable;

	void setScale(float minV, float maxV)
	{ m_scaleItrl.min = minV; m_scaleItrl.max = maxV;};
	void setRotate(float minV, float maxV)
	{ m_rotateItrl.min = minV; m_rotateItrl.max = maxV;};
	void setBias(cvS minV, cvS maxV)
	{
		m_biasItrl.resize(3);
		doF(k, 3){m_biasItrl[k].min = _f minV.val[k]; m_biasItrl[k].max = _f maxV.val[k];}
	};
	void setGain(cvS minV, cvS maxV)
	{
		m_gainItrl.resize(3);
		doF(k, 3){m_gainItrl[k].min = _f minV.val[k]; m_gainItrl[k].max = _f maxV.val[k];}
	};
};

//--------------------- Patch & PatchDistMetric ------------------------------------

class CorrSyn;

struct PatchSyn
{
	vector<cvS> pixels;
	vector<cvS> guidePixels;
};

class PatchDistMetricSyn
{
public:
	static void GetPatch(InputImageData& src, float x, float y, float s, float r, bool hr, bool vr, 
		int patchOffset, PatchSyn& patch);
	
	static float ComputePatchDist(InputImageData& src, float x, float y, CorrSyn& corr, int patchOffset);

	static float ComputePatchDistHelper(PatchSyn& vDst, PatchSyn& vSrc, cvS bias, cvS gain);
	static float CptDistDirectWithBiasAndGain(vector<cvS>& vDst, vector<cvS>& vSrc, IN cvS bias, IN cvS gain);
	static float CptDistAlphaWeight(vector<cvS>& vDst, vector<cvS>& vSrc, 
		IN cvS useAlpha, IN cvS weights, IN bool firstIsHue, IN float maxRatio, OUT cvS& alpha);
};

//--------------------- Corr & DenseCorr ---------------------------------------

class CorrSyn
{
public:
	float x, y;
	float s, r;
	float dist;
	cvS bias, gain;
	bool hr, vr; // horizenReflect, verticalReflect
	bool operator < (const CorrSyn &b) const{
		return this->dist < b.dist;
	}
	friend ostream& operator << (ostream& out, const CorrSyn& val){
		out<<"(x="<<val.x<<", y="<<val.y<<", s="<<val.s<<", r="<<val.r<<", hr = "<<val.hr<<", vr = "<<val.vr<<", dist="
			<<val.dist<<")";
		return out;
	}
	void Save(ostream& fout);
	void Load(istream& fin);
};

class DenseCorrSyn
{
public:
	DenseCorrSyn(){};
	DenseCorrSyn(int w, int h, int knn, int patchOffset = 1);
	void Identity();
	CorrSyn GetRandom(IntervalSyn& hItvl, IntervalSyn& wItvl, GPMSynRange& range, cvi* legalMask = NULL);
	void RandomInitialize(InputImageData& imgData, GPMSynRange& range);
	void UpdatePatchDistance(InputImageData& imgData);
	void ShowCorr(string imgStr);
	void ShowReflect(string imgStr);
	void ShowCorrDist(string imgStr);
	int GetCorrIdx(int r, int c){
		return (r * m_width + c) * m_knn;
	}
	CorrSyn& Get(int idx){
		return m_values[idx];
	}
	float GetDistThres(int r, int c){
		return m_values[GetCorrIdx(r, c+1) - 1].dist;
	}
	void AddCoor(int r, int c, float cx, float cy, float cs, float cr, cvS cBias, cvS cGain, bool hr, bool vr, float cdist);
	void Save(ostream& fout);
	void Load(istream& fin);
	void LevelUp(int ratio = 2);
	void AddBound(int w1, int w2, int h1, int h2);
	void HandleHoleBoundary(InputImageData& imgData, GPMSynRange& range);
public:
	int m_width, m_height, m_knn;
	int m_patchOffset;
	vector<CorrSyn> m_values; 
};

//--------------------- GPMProc ------------------------------------------------

class GPMSynProc
{

public:

	GPMSynProc(int knn = 1, int patchSize = 7, int nItr = 10, GPMSynRange* range = NULL);
	~GPMSynProc(void);

	DenseCorrSyn* RunGPM(InputImageData& src);
	void RunGPMWithInitial(InputImageData& src, DenseCorrSyn* );

	static GPMSynRange GetRange();

private:

	int Propagate(InputImageData& src, int x, int y, bool direction, 
		DenseCorrSyn& dsCor, IntervalSyn wItvl, IntervalSyn hItvl, 
		IntervalSyn swItvl, IntervalSyn shItvl, PatchSyn& dstPatch);
	int RandomSearch(InputImageData& src, int x, int y, 
		DenseCorrSyn& dsCor, IntervalSyn wItvl, IntervalSyn hItvl, PatchSyn& dstPatch);

private:

	//arguments
	int m_knn;
	int m_patchSize, m_patchOffset;
	int m_nItr;

	GPMSynRange m_range;
};

