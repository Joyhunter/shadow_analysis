#include "StdAfx.h"
#include "gpm_analysis.h"
#include "shadow_analysis.h"
#include "param_loader.h"

//---------------------------- Load Cfg --------------------------------------------------

void GPMCfg::Init()
{
	stepEnabled = true;

	//directory
	srcDir = ".//";
	resSubdir = "00_Corr//";

	//focused files
	focusedPrefix.clear();
	focusedPrefix.push_back("012");

	//Step 1: Run Grid GPM
	int patchSize_gpm = 5;
	int nItrl_gpm = 20;
	int src_colorMode_gpm = CV_BGR2Lab;
	float src_resizeRatio_gpm = 0.15f;
	float gridSizeRatio_gpm = 0.25f;
	float pymResizeRatio_gpm = 0.5f;
	float pymMinWidth_gpm = 100;
	range_gpm.setScale(0.67f,1.5f);
	range_gpm.setRotate(-1.05f * (float)CV_PI, 1.05f * (float)CV_PI);
	range_gpm.setGain(cvs(0.26, 1, 1), cvs(0.87, 1, 1));
	range_gpm.setBias(cvs(0, -6.48, -5.5), cvs(0, 8.47, 19.5));
	runGPM = true; showGPM = false;
	distThres_showGPM = 30;
}

void GPMCfg::InitFromXML(string cfgFile)
{
	tixml::XMLDoc doc;
	if(!doc.Load(cfgFile.c_str())){
		cout<<cfgFile<<" is missing. exit!\n";
		return;
	}

	stepEnabled = (doc.get<int>("gpmStep.enabled.@val", 0) == 1);

	//directory
	srcDir = doc.get<string>("srcDir.@val", 0);
	resSubdir = doc.get<string>("gpmStep.subDir.@val", 0);

	//focused files
	int num = doc.get<int>("focusedImgs.@num", 0);
	focusedPrefix.clear();
	doF(k, num)	focusedPrefix.push_back(doc.get<string>("focusedImgs.focuedPrefix.@val", k));

	//Step 1: Run Grid GPM
	patchSize_gpm = doc.get<int>("gpmStep.patchSize.@val", 0);
	nItrl_gpm = doc.get<int>("gpmStep.gpm.itrlN.@val", 0);
	src_colorMode_gpm = CV_BGR2Lab;
	src_resizeRatio_gpm = doc.get<float>("gpmStep.srcResizeRatio.@val", 0);
	gridSizeRatio_gpm = doc.get<float>("gpmStep.grid.@gridSizeRatio", 0);
	pymResizeRatio_gpm = doc.get<float>("gpmStep.pym.@resizeRatio", 0);
	pymMinWidth_gpm = doc.get<float>("gpmStep.pym.@minWidth", 0);
	string prefix = "gpmStep.gpm.";
	range_gpm.setScale(doc.get<float>(prefix + "scaleRotate.scale.@min", 0), doc.get<float>(prefix + "scaleRotate.scale.@max", 0));
	range_gpm.setRotate(doc.get<float>(prefix + "scaleRotate.rotate.@min", 0) * (float)CV_PI, 
		doc.get<float>(prefix + "scaleRotate.rotate.@max", 0) * (float)CV_PI);
	range_gpm.setGain(cvs(doc.get<float>(prefix + "gainBias.gain.@min1", 0), doc.get<float>(prefix + "gainBias.gain.@min2", 0), 
		doc.get<float>(prefix + "gainBias.gain.@min3", 0)), 
		cvs(doc.get<float>(prefix + "gainBias.gain.@max1", 0), doc.get<float>(prefix + "gainBias.gain.@max2", 0), 
		doc.get<float>(prefix + "gainBias.gain.@max3", 0)));
	range_gpm.setBias(cvs(doc.get<float>(prefix + "gainBias.bias.@min1", 0), doc.get<float>(prefix + "gainBias.bias.@min2", 0), 
		doc.get<float>(prefix + "gainBias.bias.@min3", 0)), 
		cvs(doc.get<float>(prefix + "gainBias.bias.@max1", 0), doc.get<float>(prefix + "gainBias.bias.@max2", 0), 
		doc.get<float>(prefix + "gainBias.bias.@max3", 0)));
	range_gpm.flipEnabled = (doc.get<int>(prefix + "flipEnabled.@val", 0) == 1);
	range_gpm.srEnabled = (doc.get<int>(prefix + "scaleRotate.@enabled", 0) == 1);
	if(!range_gpm.srEnabled) range_gpm.setScale(1, 1), range_gpm.setRotate(0, 0);
	range_gpm.gbEnabled = (doc.get<int>(prefix + "gainBias.@enabled", 0) == 1);
	if(!range_gpm.gbEnabled) range_gpm.setGain(cvs(1, 1, 1), cvs(1, 1, 1)), range_gpm.setBias(cvs(0, 0, 0), cvs(0, 0, 0));
	range_gpm.rsEnabled = (doc.get<int>(prefix + "randomSearch.@enabled", 0) == 1);
	range_gpm.rs_radius = doc.get<float>(prefix + "randomSearch.radius.@val", 0);
	range_gpm.rs_minRadius = doc.get<float>(prefix + "randomSearch.minRadius.@val", 0);
	range_gpm.rs_scaleFactor = doc.get<float>(prefix + "randomSearch.scaleFactor.@val", 0);
	range_gpm.nCores = doc.get<int>("gpmStep.coreN.@val", 0);
	runGPM = (doc.get<int>("gpmStep.runGPM.@val", 0) == 1);
	showGPM = (doc.get<int>("gpmStep.showGPM.@val", 0) == 1);
	distThres_showGPM = doc.get<float>("gpmStep.showGPM.distThres.@val", 0);
}

void VoteCfg::Init()
{
	stepEnabled = true;

	//directory
	srcDir = ".//";
	corrDir = "00_Corr//";
	resSubdir = "01_PredictParam//";

	//focused files
	focusedPrefix.clear();
	focusedPrefix.push_back("012");

	useNaive = true;

	patchSize = 5;
	distThres = 20;
	useTopN = 0;

	isTraining = false;
	modelFile = "..//rt.xml";
	rtParam.matchNPerPatch = 12;
	rtParam.sampleStep = 1;
	rtParam.max_depth = 25;
	rtParam.min_sample_count = 5;
	rtParam.reggression_accurancy = 0;
	rtParam.max_categories = 15;
	rtParam.nactive_vars = 4;
	rtParam.max_tree_N = 100;
	rtParam.forest_accurancy = 0.01f;
	rtParam.forest_accurancyAB = 1;
	gtSubDir = "00_GT//";

	usePrediction = false;
}

void VoteCfg::InitFromXML(string cfgFile)
{
	tixml::XMLDoc doc;
	if(!doc.Load(cfgFile.c_str())){
		cout<<cfgFile<<" is missing. exit!\n";
		return;
	}

	stepEnabled = (doc.get<int>("voteStep.enabled.@val", 0) == 1);

	//directory
	srcDir = doc.get<string>("srcDir.@val", 0);
	corrDir = doc.get<string>("gpmStep.subDir.@val", 0);
	resSubdir = doc.get<string>("voteStep.subDir.@val", 0);

	//focused files
	int num = doc.get<int>("focusedImgs.@num", 0);
	focusedPrefix.clear();
	doF(k, num)	focusedPrefix.push_back(doc.get<string>("focusedImgs.focuedPrefix.@val", k));

	useNaive = (doc.get<int>("voteStep.useNaive.@val", 0) == 1);

	//naive voting
	patchSize = doc.get<int>("voteStep.naive.patchSize.@val", 0);
	distThres = doc.get<float>("voteStep.naive.distThres.@val", 0);
	useTopN = doc.get<int>("voteStep.naive.useTopN.@val", 0);

	//prediction
	isTraining = (doc.get<int>("voteStep.training.enabled.@val", 0) == 1);
	modelFile = doc.get<string>("voteStep.training.model.@fileDir", 0);
	string prefix = "voteStep.training.rt.";
	rtParam.matchNPerPatch = doc.get<int>(prefix + "matchNPerPatch.@val", 0);
	rtParam.sampleStep = doc.get<int>(prefix + "sampleStep.@val", 0);
	rtParam.max_depth = doc.get<int>(prefix + "max_depth.@val", 0);
	rtParam.min_sample_count = doc.get<int>(prefix + "min_sample_count.@val", 0);
	rtParam.reggression_accurancy = doc.get<float>(prefix + "reggression_accurancy.@val", 0);
	rtParam.max_categories = doc.get<int>(prefix + "max_categories.@val", 0);
	rtParam.nactive_vars = doc.get<int>(prefix + "nactive_vars.@val", 0);
	rtParam.max_tree_N = doc.get<int>(prefix + "max_tree_N.@val", 0);
	rtParam.forest_accurancy = doc.get<float>(prefix + "forest_accurancy.@val", 0);
	rtParam.forest_accurancyAB = doc.get<float>(prefix + "forest_accurancyAB.@val", 0);
	gtSubDir = doc.get<string>("gtStep.resSubDir.@val", 0);

	usePrediction = (doc.get<int>("voteStep.usePrediction.@val", 0) == 1);
}

//----------------------------- Init functions --------------------------------------------

GPMAnalysisProc::GPMAnalysisProc(void)
{
}

GPMAnalysisProc::~GPMAnalysisProc(void)
{
}

GPMCfg GPMAnalysisProc::cfg;
VoteCfg GPMAnalysisProc::vcfg;

void GPMAnalysisProc::LoadCfg(string cfgFile)
{
	cfg.InitFromXML(cfgFile);
	vcfg.InitFromXML(cfgFile);
}

//----------------------------- interface 1: Run Grid GPM for all images -------------------

void GPMAnalysisProc::RunGPMForAllImages()
{
	COutYel("-- Grid GPM --\n");
	wMkDir(cfg.srcDir + cfg.resSubdir);
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

	doFv(i, m_imgNames)
	{
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		if(cfg.focusedPrefix.size() > 0 && find(cfg.focusedPrefix.begin(), cfg.focusedPrefix.end(), imgPrefix) == cfg.focusedPrefix.end())
			continue;
		COutTeal("Processing " + imgPrefix + ".png..\n");

		ImgContainer img(cfg.srcDir + fileName, cfg.src_resizeRatio_gpm, cfg.src_colorMode_gpm);

		string corrFileName = cfg.srcDir + cfg.resSubdir + imgPrefix + ".corr";

		LmnIvrtPatchDistMetric metric2;
		GridGPMProc proc(&metric2, cfg.gridSizeRatio_gpm, 1, cfg.patchSize_gpm, cfg.nItrl_gpm, &(cfg.range_gpm));

		if(cfg.runGPM) proc.RunGridGPMMultiScale(img, corrFileName, cfg.pymResizeRatio_gpm, cfg.pymMinWidth_gpm);
		//if(img.srcR() == NULL) img.GenerateResizedImg(1);
		if(cfg.showGPM) proc.ShowGPMResUI(img, corrFileName, cfg.distThres_showGPM);
	}
}

//-------------------------- interface 2: naive vote ------------------------------------

void GPMAnalysisProc::RunNaiveVoting()
{
	COutYel("-- Naive Voting --\n");
	wMkDir(vcfg.srcDir + vcfg.resSubdir);
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

// 	omp_set_num_threads(6);
// #pragma omp parallel for
	doF(k, _i m_imgNames.size())
	{
		string fileName = m_imgNames[k];
		string imgPrefix = m_imgNames[k].substr(0, m_imgNames[k].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		if(vcfg.focusedPrefix.size() > 0 && find(vcfg.focusedPrefix.begin(), vcfg.focusedPrefix.end(), imgPrefix) == vcfg.focusedPrefix.end())
			continue;
		COutTeal("Processing " + imgPrefix + ".png..\n");

		//read corr file
		string corrFile = vcfg.srcDir + vcfg.corrDir + imgPrefix + ".corr";
		DenseCorrBox2D box;
		box.Load(corrFile);
		if(box.m_dstH == 0)
		{
			COutRed("Coor file " + corrFile + " not exist error.\n");
			continue;
		}
		int w = box.m_dstW, h = box.m_dstH;

		//init image
		ImgContainer img(vcfg.srcDir + fileName, 1.0, cfg.src_colorMode_gpm);
		img.GenerateResizedImg(w, h);

		cvi* param = cvci323(w, h);
		doFcvi(img.srcR(), i, j)
		{
			MultiCorr corrs = box.GetCorrsPerGrid(i, j);

			vector<double> gbV(6);

			Patch srcPatch, dstPatch;
			PatchDistMetric::GetPatch(img, _f i, _f j, 1.f, 0.f, false, false, (vcfg.patchSize-1)/2, srcPatch);
			
			vector<vector<double> > gbVs(0);
			vector<pair<double, int> > Ls(0); 

			doFv(k, corrs)
			{
				Corr& v = corrs[k];
				if(v.dist > vcfg.distThres) continue;
				PatchDistMetric::GetPatch(img, v.x, v.y, v.s, v.r, v.hr, v.vr, (vcfg.patchSize-1)/2, dstPatch);
				double dist;
				double L = OptimizeGbVLab(srcPatch, dstPatch, gbV, dist);
				Ls.push_back(make_pair(L, Ls.size()));
				gbVs.push_back(gbV);
			}

			if(Ls.size() == 0)
			{
				doF(k, 3) gbV[2*k] = 1;
				doF(k, 3) gbV[2*k+1] = 0;
			}
			else
			{
				sort(Ls.begin(), Ls.end(), [](const pair<double, int>& v1, const pair<double, int>& v2){
					return v1.first > v2.first;
				});
				int idx = min2(vcfg.useTopN, _i Ls.size() - 1);
				gbV = gbVs[Ls[idx].second];
			}

			cvs2(param, i, j, cvs(gbV[0], gbV[3], gbV[5]));
		}

		ParamLoader::SaveParamToDir(vcfg.srcDir + vcfg.resSubdir, imgPrefix, param);
		ParamLoader::ShowParamInDir(vcfg.srcDir + vcfg.resSubdir, imgPrefix, param);
		cvri(param);
	}
}

double GPMAnalysisProc::OptimizeGbVLab(Patch& srcPatch, Patch& dstPatch, vector<double>& gbV, double& dist)
{
	cvS srcAvg = cvs(0, 0, 0), dstAvg = cvs(0, 0, 0);
	doFv(k, srcPatch.hls)
	{
		srcAvg += srcPatch.hls[k];
	}
	srcAvg /= srcPatch.hls.size();
	doFv(k, dstPatch.hls)
	{
		dstAvg += dstPatch.hls[k];
	}
	dstAvg /= dstPatch.hls.size();

	gbV[0] = srcAvg.val[0] / dstAvg.val[0];
	gbV[1] = 0;
	gbV[2] = 1;
	gbV[3] = - srcAvg.val[1] + dstAvg.val[1];
	gbV[4] = 1;
	gbV[5] = - srcAvg.val[2] + dstAvg.val[2];

	cvS gain = cvs(gbV[0], gbV[2], gbV[4]);
	cvS bias = cvs(gbV[1], gbV[3], gbV[5]);
	float sum = 0;
	doFv(i, srcPatch.hls)
	{
		sum += _f cvSDSqr(srcPatch.hls[i], dstPatch.hls[i]*gain+bias);
	}
	dist = sqrt(sum / srcPatch.hls.size());

	return dstAvg.val[0];
}

//---------------------------------------------------------------------------------------

//5: generate GT 

//todo 1: compute matching accurancy
//2: vote: / prediction
//3: voteError
//4: train model

void GPMAnalysisProc::SetFileDir(string fileDir)
{
	//cfg.srcDir = fileDir;
}

void GPMAnalysisProc::ReadParam(string cfgFile)
{
// 	ifstream fin(cfgFile.c_str());
// 	float v1, v2, v3, v4, v5, v6;
// 	fin>>v1>>v2;
// 	cfg.range_gpm.setScale(v1, v2);
// 	fin>>v1>>v2;
// 	cfg.range_gpm.setRotate(v1 * (float)CV_PI, v2 * (float)CV_PI);
// 	fin>>v1>>v2>>v3>>v4>>v5>>v6;
// 	cfg.range_gpm.setGain(cvs(v1, v2, v3), cvs(v4, v5, v6)); //0.92
// 	fin>>v1>>v2>>v3>>v4>>v5>>v6;
// 	cfg.range_gpm.setBias(cvs(v1, v2, v3), cvs(v4, v5, v6)); //21.4
// 	fin>>param.downRatio;
// 	fin>>param.multiLevelRatio;
// 	fin>>param.extraLevelN;
// 	fin>>param.patchSize;
// 	fin>>param.distThres;
// 	fin>>rangeDir;
// 	fin.close();
}

string GPMAnalysisProc::GetCorrFileDir(string imgName)
{
// 	return "corr//" + rangeDir + imgName.substr(0, imgName.size()-4)
// 		+ "_" + toStr(param.downRatio) + "_" + toStr(param.multiLevelRatio) + "_" + toStr(param.extraLevelN)
// 		+ "_" + toStr(param.patchSize) + ".corr";
	return "";
}







void GPMAnalysisProc::ShdwAnlysis()
{
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

	//0.01: (0.29, 0.94). (-5.47, 8.47). (-6.47461, 21.416)
	//0.05: (0.391602, 0.827148). (-3.48633, 6.47461). (-3.48633, 18.4277).
	//0.1:  (0.43457, 0.780273). (-2.49023, 4.48242). (-2.49023, 17.4316).
	//0.15: (0.463867, 0.743164). (-1.49414, 4.48242). (-0.498047, 16.4355).
	//0.2: (0.493164, 0.717773). (-1.49414, 3.48633). (0.498047, 15.4395).
	//0.25: (0.510742, 0.698242). (-0.498047, 3.48633). (4.48242, 14.4434).
	//0.3: (0.522461, 0.682617). (-0.498047, 2.49023). (6.47461, 13.4473).
	//0.4: (0.549805, 0.649414). (0.498047, 2.49023). (7.4707, 11.4551).
	//0.49: (0.594727, 0.608398). (1.49414, 1.49414). (9.46289, 9.46289).

	//range 1(0.01): 0.26 0.87  -6.48 8.47   -5.5 19.5
	//range 2(0.1): 0.43 0.78  -2.49 4.48   -2.49 17.43
	//range 3(0.15): 0.46 0.74  -1.49 4.49   -0.49 16.44
	//range 4(0.2): 0.49 0.72  -1.49 3.49  0.49 15.44  set range wrong..
	//range 5(0.25): 0.51 0.70  -0.49 3.49  4.48 14.44
	//range 6(0.3): 0.52 0.68  -0.49 2.49  6.47 13.45
	//range 7(0.4): 0.54 0.65  0.49 2.49  7.47 11.46
	//range 8(0.49): 0.59 0.61  1.49 1.50  9.46 9.47 
	//range 9(0.5): 0.60 0.60  1.50 1.50  9.46 9.46

	cfg.range_gpm.setScale(0.67f, 1.5f);
	cfg.range_gpm.setRotate(-1.05f * (float)CV_PI, 1.05f * (float)CV_PI);
	cfg.range_gpm.setGain(cvs(0.6, 1, 1), cvs(0.6, 1, 1)); //0.92
	cfg.range_gpm.setBias(cvs(0, 1.5, 9.46), cvs(0, 1.5, 9.46)); //21.4
// 	param.downRatio = 6; param.multiLevelRatio = 2; param.extraLevelN = 0;
// 	param.patchSize = 7;
// 	param.colorMode = CV_BGR2Lab;
// 	param.distThres = 30;

// 	rangeDir = "range4//";
// 	parVoteDir = "VoteParam//";
// 	parGTDir = "GTParam6//";

	//ReadParam();

// 	if(1)
// 	{
// 		//RunGPMForAllImages(); return;
// 
// 		ofstream fout("sta.txt");
// 		doF(i, 128)
// 		{
// 			cout<<"\r"<<i;
// 			param.distThres = _f 2 * i;
// 			GetStatistics(fout);
// 		}
// 		fout.close();
// 	}

	test(); return;

// 	VoteInitMask(); return;
// 	ofstream fout2("voteError.txt");
// 	doF(i, 40)
// 	{
// 		cout<<"\r"<<i*2;
// 		param.distThres = _f 2 * i;
// 		VoteInitMask(); 
// 		fout2<<GetVoteError()<<endl;
// 	}
// 	fout2.close();
// 	return;
}

void GPMAnalysisProc::test()
{
	ifstream fin("vote.cfg");
// 	fin>>rangeDir;
// 	fin>>parVoteDir;
// 	fin>>parGTDir;
// 	fin>>param.downRatio;
// 	fin>>param.multiLevelRatio;
// 	fin>>param.extraLevelN;
// 	fin>>param.distThres;
	fin.close();
	wMkDir(cfg.srcDir + parVoteDir);
	//VoteInitMask(); GetVoteError(); pause; return;
	ofstream fout2("voteError.txt");
	doF(i, 40)
	{
		//if(2*i < param.distThres) continue;
		cout<<"\r"<<i*2;
		//param.distThres = _f 2 * i;
		RunNaiveVoting(); 
		cvS res = GetVoteError();
		fout2<<res.val[0]<<" "<<res.val[1]<<" "<<res.val[2]<<endl;
	}
	fout2.close();

}

string GPMAnalysisProc::GetShadowSegDir(string imgName)
{
	return "mask//" + imgName.substr(0, imgName.size()-4) + "s.png";
}
string GPMAnalysisProc::GetMaterialSegDir(string imgName)
{
	return "mask//" + imgName.substr(0, imgName.size()-4) + "m.png";
}


void GPMAnalysisProc::GetStatistics(ofstream& fout)
{
// 	GPMStatis sta;
// 
// 	//cout<<"\nBegin analysis matching result for all images...\n";
// 	omp_set_num_threads(6);
// #pragma omp parallel for
// 	doF(i, _i m_imgNames.size())
// 	{
// 		if(m_imgNames[i][m_imgNames[i].size()-5] == 'n') continue;
// 		//cout<<"--Analysis image "<<m_imgNames[i]<<"...\n";
// 		GPMStatis sta1;
// 
// 		string corrFile = GetCorrFileDir(m_imgNames[i]);
// 		string shadowMaskFile = GetShadowSegDir(m_imgNames[i]);
// 		string matMaskDir = GetMaterialSegDir(m_imgNames[i]);
// 
// 		cvi* shMaskOri = cvlig(cfg.srcDir + shadowMaskFile);
// 		cvi* shMask = cvci81(shMaskOri->width / param.downRatio, shMaskOri->height / param.downRatio);
// 		cvResize(shMaskOri, shMask, CV_INTER_NN);
// 		//cvsi(shMask);
// 		cvi* matMaskOri = cvlig(cfg.srcDir + matMaskDir);
// 		cvi* matMask = cvci(shMask);
// 		cvResize(matMaskOri, matMask, CV_INTER_NN);
// 		cvri(shMaskOri); cvri(matMaskOri);
// 
// 		DenseCorrBox2D box;
// 		box.Load(cfg.srcDir + corrFile);
// 
// 		doFcvi(shMask, i, j)
// 		{
// 			MultiCorr corrs = box.GetCorrsPerGrid(i, j);
// 
// 			int shdV = _i cvg20(shMask, i, j), matV = _i cvg20(matMask, i, j);
// 			
// 			//unshadow
// 			if(shdV > 0)
// 			{
// 				doFv(k, corrs)
// 				{
// 					Corr& cor = corrs[k];
// 					sta1.unShdwMatchN ++;
// 					if(cor.dist <= param.distThres)
// 					{
// 						sta1.wum++;
// 					}
// 				}
// 				continue;
// 			}
// 
// 			bool firstMatch = false;
// 			sta1.pixelN++;
// 			doFv(k, corrs)
// 			{
// 				Corr& cor = corrs[k];
// 				sta1.matchN++;
// 				if(cor.dist > param.distThres)
// 				{
// 					sta1.wm++;
// 					continue;
// 				}
// 				int shdV2 = _i cvg20(shMask, round(cor.x), round(cor.y));
// 				int matV2 = _i cvg20(matMask, round(cor.x), round(cor.y));
// 				if(shdV != shdV2 && matV == matV2)
// 				{
// 					if(firstMatch == false)
// 					{
// 						sta1.rp++;
// 						firstMatch = true;
// 					}
// 					sta1.dssm++;
// 				}
// 				else if(shdV != shdV2 && matV != matV2) sta1.dsdm++;
// 				else if(shdV == shdV2 && matV == matV2){
// 					//cout<<i<<" "<<j<<" "<<cor.x<<" "<<cor.y<<" "<<shdV<<" "<<shdV2<<" "<<cor.dist<<" "<<endl; pause;
// 					sta1.sssm++;
// 				}
// 				else sta1.ssdm++;
// 			}
// 		}
// 		//sta1.print();
// 		sta.add(sta1);
// 
// 		cvri(shMask); cvri(matMask);
// 	}
// 	//cout<<"--Analysis matching result for all images complete! Summary: "<<endl;
// 	fout<<param.distThres<<"\t";
// 	sta.print2(fout);
// 	fout<<endl;
}

cvS GPMAnalysisProc::GetVoteError()
{
	cout<<"Begin calculating voting error using GTParam and VoteParam...\n";
	cvS allError = cvs(0, 0, 0);
	doFv(i, m_imgNames)
	{
		if(m_imgNames[i][m_imgNames[i].size()-5] == 'n') continue;
		cout<<"\r--Analysis image "<<m_imgNames[i]<<"...";
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 4);

		string voteResDirLGain = cfg.srcDir + parVoteDir + imgPrefix + "_0_gain.txt";
		string gtResDirLGain = cfg.srcDir + parGTDir + imgPrefix + "_0_gain.txt";
		string voteResDirABias = cfg.srcDir + parVoteDir + imgPrefix + "_1_bais.txt";
		string gtResDirABias = cfg.srcDir + parGTDir + imgPrefix + "_1_bais.txt";
		string voteResDirBBias = cfg.srcDir + parVoteDir + imgPrefix + "_2_bais.txt";
		string gtResDirBBias = cfg.srcDir + parGTDir + imgPrefix + "_2_bais.txt";

		cvS sumError = cvs(0, 0, 0);
		ifstream fin(voteResDirLGain.c_str()), fin2(gtResDirLGain.c_str());
		ifstream fin3(voteResDirABias.c_str()), fin4(gtResDirABias.c_str());
		ifstream fin5(voteResDirBBias.c_str()), fin6(gtResDirBBias.c_str());
		int w, h, w2, h2;
		fin3>>w>>h; fin4>>w>>h; fin5>>w>>h; fin6>>w>>h; fin>>w>>h; fin2>>w2>>h2; 
		if(w != w2 || h != h2)
		{
			cout<<"Error!"<<endl;
		}
		else
		{
			double v1, v2, v3, v4, v5, v6;
			doF(k, w*h)
			{
				fin>>v1; fin2>>v2; fin3>>v3; fin4>>v4; fin5>>v5; fin6>>v6;
				sumError += cvs(abs(v1 - v2), abs(v3-v4), abs(v5-v6));
			}
			sumError /= w*h;
		}
		fin.close(); fin2.close();
		allError += sumError;
	}
	allError /= m_imgNames.size() / 2;
	cout<<"\rVote Error computing complete, vote error = "<<allError.val[0]<<" "<<allError.val[1]<<" "
		<<allError.val[2]<<".\n";
	return allError;
}
