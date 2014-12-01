#pragma once
#include "../shadow_removal/gpm_proc.h"

struct GPMParam
{
	GPMRange range;

	int downRatio, multiLevelRatio, extraLevelN; 
	int patchSize;

	int colorMode;

	float distThres;
};

struct RTparam
{
	int matchNPerPatch;
	int sampleStep;
	int max_depth;
	int min_sample_count;
	float reggression_accurancy;
	int max_categories;
	int nactive_vars;
	int max_tree_N;
	float forest_accurancy;
	float forest_accurancyAB;
};

struct GPMStatis
{
	int matchN; //all matches
	int dssm; //different shadow, same material
	int sssm; //same shadow, same material
	int dsdm; //different shadow, different material
	int ssdm; //same shadow, different material
	int wm; // wrong match

	int pixelN; // all pixels
	int rp; // right matched pixels

	int unShdwMatchN; // all unshadow matches
	int wum; // wrong unshadow matches

	GPMStatis():matchN(0), dssm(0), sssm(0), dsdm(0), ssdm(0), wm(0), pixelN(0), rp(0), unShdwMatchN(0), wum(0){};
	void add(GPMStatis& v)
	{
		matchN += v.matchN; dssm += v.dssm; sssm += v.sssm; dsdm += v.dsdm; ssdm += v.ssdm; wm += v.wm;
		pixelN += v.pixelN; rp += v.rp;
		unShdwMatchN += v.unShdwMatchN; wum += v.wum;
	}
	void print()
	{
		cout<<"MatchN = "<<matchN<<".\n";
		cout<<"Un shadow, same material\t= "<<dssm<<"("<<((matchN == 0)?0:(_f dssm*100/matchN))<<"%).\n";
		cout<<"In shadow, same material\t= "<<sssm<<"("<<((matchN == 0)?0:(_f sssm*100/matchN))<<"%).\n";
		cout<<"Un shadow, dif material\t\t= "<<dsdm<<"("<<((matchN == 0)?0:(_f dsdm*100/matchN))<<"%).\n";
		cout<<"In shadow, dif material\t\t= "<<ssdm<<"("<<((matchN == 0)?0:(_f ssdm*100/matchN))<<"%).\n";
		cout<<"Wrong match\t\t\t= "<<wm<<"("<<((matchN == 0)?0:(_f wm*100/matchN))<<"%).\n";
		cout<<"PixelN = "<<pixelN<<".\n";
		cout<<"Right matched pixel\t\t= "<<rp<<"("<<((pixelN == 0)?0:(_f rp*100/pixelN))<<"%).\n";
		cout<<"UnShdwMatchN = "<<unShdwMatchN<<".\n";
		cout<<"Wrong unshadow matchs\t\t= "<<wum<<"("<<((unShdwMatchN == 0)?0:(_f wum*100/unShdwMatchN))<<"%).\n";
	}
	void print2(ofstream& fout)
	{
 		fout<<((matchN == 0)?0:(_f dssm*100/matchN))<<"\t";
 		fout<<((pixelN == 0)?0:(_f rp*100/pixelN))<<"\t";
 		fout<<((unShdwMatchN == 0)?0:(_f wum*100/unShdwMatchN))<<"\t";
		fout<<(((matchN - wm + wum) == 0)?100:(_f dssm*100 / (matchN - wm + wum)));
	}
};

//----------------------- Cfg of Grid-GPM proc --------------------------
struct GPMCfg
{
	//enabled
	bool stepEnabled;

	//directory
	string srcDir;
	string resSubdir;

	//focused files
	vector<string> focusedPrefix;

	//Step 1: Run Grid GPM
	int patchSize_gpm;
	int nItrl_gpm;
	int src_colorMode_gpm;
	float src_resizeRatio_gpm;
	float gridSizeRatio_gpm;
	float pymResizeRatio_gpm;
	float pymMinWidth_gpm;
	GPMRange range_gpm;
	bool runGPM, showGPM;
	float distThres_showGPM;

	void Init();
	void InitFromXML(string cfgFile);
};

struct VoteCfg
{
	bool stepEnabled;

	string srcDir;
	string corrDir;
	string resSubdir;

	//focused files
	vector<string> focusedPrefix;

	//naive vote
	bool useNaive;
	int patchSize;
	float distThres;
	int useTopN;

	//training
	bool isTraining;
	string modelFile;
	RTparam rtParam;
	string gtSubDir;

	bool usePrediction;

	void Init();
	void InitFromXML(string cfgFile);
};

//------------------------------------------------------------------------

class GPMAnalysisProc
{

public:
	
	GPMAnalysisProc(void);
	~GPMAnalysisProc(void);

	void LoadCfg(string cfgFile);

	//Run gpm interface: Run Grid GPM for all images, and save file to corr
	void RunGPMForAllImages();

	//Naive voting
	void RunNaiveVoting();

	//Predict using r trees
	void PredictUseRTrees();
	
	//vote learning
	void TrainVoteRTrees();

	static GPMCfg cfg;
	static VoteCfg vcfg;

private:

	//Get corr file save dir of imgName
	string GetCorrFileDir(string imgName);


	//abandoned
	void InitRTParam();
public:

	void SetFileDir(string fileDir);

	void ShdwAnlysis();

	//learn to predict
	void VoteAnalysis();

	static cvi* LoadGTFile(string gtDir, string imgPrefix);

private:

	//Read cfg to GPMParam param
	void ReadParam(string cfgFile = "gpm.cfg");

	void GetStatistics(ofstream& fout);

	double OptimizeGbVLab(Patch& srcPatch, Patch& dstPatch, vector<double>& gbV, double& dist);
	cvS GetVoteError();

	void test();

	string GetShadowSegDir(string imgName);
	string GetMaterialSegDir(string imgName);

private:

	//string m_fileDir;
	vector<string> m_imgNames;

	//param
	//GPMParam param;

	//string rangeDir;
	string parVoteDir;
	string parGTDir;

	//string parPredictDir;
};

