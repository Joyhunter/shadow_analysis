#include "StdAfx.h"
#include "detection_proc.h"
#include "param_loader.h"
#include "../shadow_removal/mrf_proc.h"
#include "recovery_proc.h"
#include "local_propagation/local_editprop_proc.h"
#include "alpha_matting/alpha_matting_proc.h"

//----------------------------- Detection Configure ----------------------------

void DetectionCfg::Init()
{
	stepEnabled = true;

	//directory
	srcDir = ".//";
	paramSubdir = "01_VoteParam//";
	resSubdir = "02_VoteDetection//";

	//focused files
	focusedPrefix.clear();
	focusedPrefix.push_back("012");

	//step 1: mrf smooth
	enabled_mrf = true;
	nLabels_mrf = 16;
	nCh_mrf = 3;
	colorsigma_mrfSTerm = 0.2f;
	sWeight_mrfSTerm = 2.0f;

	//step 2: ana
	dCfg.hishN = 64;
	dCfg.histMinV = 0.1f;
	dCfg.adjacentRatio = 1.0f;
	dCfg.shadowDegreeThres = 200.f;
	dCfg.pixelNumThres = 0.2f;
	dCfg.boxRatioThres = 0.2f;
	dCfg.cellAdjacentRadiusRatio = 0.02f;
	dCfg.regionPredictor = false;
	dCfg.regionPredictorDir = "..//";
	dCfg.regionAlphaEstimateHigherVal = 0.2f;
	dCfg.boundaryWidth = -1;
	dCfg.nLabels = 16;
	guidInitMaxRatio_data = 1.2f;
	ratioSigma_data = 0.05f; 
	whiteWeight_data = 10.0f;
	colorSigma_smooth = 0.2f;
	weight_smooth = 10.0f;

	//step 3: matting
	shdwDegreeThres_matting = 0.8f;
	maskErodeTimes_matting = 3;
	maskDilateTimes_matting = 5;

	//debug
	debugImgsOutput = true;
	debugImgsOutputDir = "result//";
	dCfg.debugImgsOutput = true;
	dCfg.debugImgsOutputDir = "result//";
}

void DetectionCfg::InitFromXML(string cfgFile)
{
	tixml::XMLDoc doc;
	if(!doc.Load(cfgFile.c_str())){
		cout<<cfgFile<<" is missing. exit!\n";
		return;
	}

	stepEnabled = (doc.get<int>("detectionStep.enabled.@val", 0) == 1);

	//directory
	srcDir = doc.get<string>("srcDir.@val", 0);
	paramSubdir = doc.get<string>("voteStep.subDir.@val", 0); 
	resSubdir = doc.get<string>("detectionStep.subDir.@val", 0);

	//focused files
	int num = doc.get<int>("focusedImgs.@num", 0);
	focusedPrefix.clear();
	doF(k, num)	focusedPrefix.push_back(doc.get<string>("focusedImgs.focuedPrefix.@val", k));

	//step 1: mrf smooth
	enabled_mrf = (doc.get<int>("detectionStep.mrfSmooth.enabled.@val", 0) == 1);
	nLabels_mrf = doc.get<int>("detectionStep.mrfSmooth.nLabels.@val", 0);
	nCh_mrf = doc.get<int>("detectionStep.mrfSmooth.nChannels.@val", 0);
	colorsigma_mrfSTerm = doc.get<float>("detectionStep.mrfSmooth.smoothTerm.colorSigma.@val", 0);
	sWeight_mrfSTerm = doc.get<float>("detectionStep.mrfSmooth.smoothTerm.weight.@val", 0);

	//step 2: ana
	dCfg.hishN = doc.get<int>("detectionStep.analysis.hist.histN.@val", 0);
	dCfg.histMinV = doc.get<float>("detectionStep.analysis.hist.histMinV.@val", 0);
	dCfg.adjacentRatio = doc.get<float>("detectionStep.analysis.hist.adjacentRatio.@val", 0);
	dCfg.shadowDegreeThres = doc.get<float>("detectionStep.analysis.cellFilter.shadowDegreeThres.@val", 0);
	dCfg.pixelNumThres = doc.get<float>("detectionStep.analysis.cellFilter.pixelNumThres.@val", 0);
	dCfg.boxRatioThres = doc.get<float>("detectionStep.analysis.cellFilter.boxRatioThres.@val", 0);
	dCfg.cellAdjacentRadiusRatio = doc.get<float>("detectionStep.analysis.cellAdjacentRadiusRatio.@val", 0);
	dCfg.regionPredictor = (doc.get<int>("detectionStep.analysis.regionPredictor.@enabled", 0) == 1);
	dCfg.regionPredictorDir = doc.get<int>("detectionStep.analysis.regionPredictor.@val", 0);
	dCfg.regionAlphaEstimateHigherVal = doc.get<float>("detectionStep.analysis.mrf.regionAlphaEstimate.@higherVal", 0);
	dCfg.boundaryWidth = doc.get<int>("detectionStep.analysis.mrf.boundaryWidth.@val", 0);
	dCfg.nLabels = doc.get<int>("detectionStep.analysis.mrf.nLabels.@val", 0);
	guidInitMaxRatio_data = doc.get<float>("detectionStep.analysis.mrf.dataTerm.@guidInitMaxRatio", 0);
	ratioSigma_data = doc.get<float>("detectionStep.analysis.mrf.dataTerm.@ratioSigma", 0);
	whiteWeight_data = doc.get<float>("detectionStep.analysis.mrf.dataTerm.@whiteWeight", 0);
	colorSigma_smooth = doc.get<float>("detectionStep.analysis.mrf.smoothTerm.@colorSigma", 0);
	weight_smooth = doc.get<float>("detectionStep.analysis.mrf.smoothTerm.@weight", 0);

	//step 3: matting
	shdwDegreeThres_matting = doc.get<float>("detectionStep.matting.shdwDegreeThres.@val", 0);
	maskErodeTimes_matting = doc.get<int>("detectionStep.matting.maskErodeTimes.@val", 0);
	maskDilateTimes_matting = doc.get<int>("detectionStep.matting.maskDilateTimes.@val", 0);

	//debug
	debugImgsOutput = (doc.get<int>("detectionStep.debug.debugResult.@output", 0) == 1);
	debugImgsOutputDir = doc.get<string>("detectionStep.debug.debugResult.@dir", 0);
	dCfg.debugImgsOutput = debugImgsOutput;
	dCfg.debugImgsOutputDir = debugImgsOutputDir;

	// 	cout<<patchsize_naiveRecovery<<" "<<shdwDegreeThres_maskGenerate<<" "<<shdwBoundDilateTimes_maskGenerate<<" "<<patchRadius_localColorCorrection
	// 		<<" "<<correctionStepRatio_localColorCorrection<<" "<<pyramidExtraLevels_pyramid<<" "<<pyramidResizeRatio_pyramid<<" "
	// 		<<patchRadius_old;
}

DetectionCfg DetectionProc::cfg;

DetectionProc::DetectionProc(void)
{
}

DetectionProc::~DetectionProc(void)
{
}

void DetectionProc::LoadCfg(string cfgFile)
{
	cfg.InitFromXML(cfgFile);
}

//---------------------------------------------------------------------------------------

void DetectionProc::DetectCastShadow()
{
	COutYel("-- Detection --\n");
	wMkDir(cfg.srcDir + cfg.resSubdir);
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

	doFv(i, m_imgNames)
	{
		//get image prefix
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		if(cfg.focusedPrefix.size() > 0 && find(cfg.focusedPrefix.begin(), cfg.focusedPrefix.end(), imgPrefix) == cfg.focusedPrefix.end())
			continue;
		COutTeal("Processing " + imgPrefix + ".png..\n");

		//load source image and param prediction
		cvi* srcImg = cvlic(cfg.srcDir + imgPrefix + ".png");
		cvi* param = ParamLoader::LoadParamFromDir(cfg.srcDir + cfg.paramSubdir, imgPrefix);//m_srcDir

		//cast shadow detect
		cvi* paramNew = DetectCastShadow(srcImg, param);

		//save detect result
		ParamLoader::SaveParamToDir(cfg.srcDir + cfg.resSubdir, imgPrefix, paramNew);
		ParamLoader::ShowParamInDir(cfg.srcDir + cfg.resSubdir, imgPrefix, paramNew);

		cvri(srcImg); cvri(param); cvri(paramNew);
	}
}

cvi* DetectionProc::DetectCastShadow(cvi* srcImg, cvi* param)
{
	cvi* paramResize = cvci323(srcImg);
	cvResize(param, paramResize);
	if(cfg.debugImgsOutput)
	{
		ParamLoader::ShowParamInDir(cfg.debugImgsOutputDir, "step0_init", paramResize);
		cvi* resImg = RecoverProc::ImgRecoverNaive(srcImg, paramResize, 1);
		cvsi(cfg.debugImgsOutputDir + "res_step0.png", resImg); cvri(resImg);
	}

	COutput("Start step 1: overall MRF Smooth..\n", CC_DARKGREEN);
	//step 1: Up-sampling, overall smooth
	//srcImg = cvlic("result//004blur2.png");
	//cvi* smoothParam = cvci(param);
	cvi* smoothParam = paramResize;
	if(cfg.enabled_mrf) smoothParam = MRFSmooth(srcImg, paramResize, cfg.nLabels_mrf);
	if(cfg.debugImgsOutput)
	{
		ParamLoader::ShowParamInDir(cfg.debugImgsOutputDir, "step1_mrfSmooth", smoothParam);
		cvi* resImg = RecoverProc::ImgRecoverNaive(srcImg, smoothParam, 1);
		cvsi(cfg.debugImgsOutputDir + "res_step1.png", resImg); cvri(resImg);
	}
	COutGreen("\r                                   \rStep 1: MRF smooth res got.\n");
	
	COutput("Start step 2: shadow decompose, cast shadow detection..\n", CC_DARKGREEN);
	//step 2: Analysis, decompose, another MRF smooth
	smoothParam = DecmpsSmooth(srcImg, smoothParam, 16);
	if(cfg.debugImgsOutput)
	{
		ParamLoader::ShowParamInDir(cfg.debugImgsOutputDir, "step2_decpmSmooth", smoothParam);
		cvi* resImg = RecoverProc::ImgRecoverNaive(srcImg, smoothParam, 1);
		cvsi(cfg.debugImgsOutputDir + "res_step2.png", resImg); cvri(resImg);
	}
	COutGreen("\r                                   \rStep 2: Cast shadow got.\n");

	COutput("Start step 3: matting shadow boundary..", CC_DARKGREEN);
	//step 3: Matting(trick)
	//cvi* smoothSrc = cvci(srcImg);
	//cvSmooth(srcImg, smoothSrc, 2, 31, 31);
	smoothParam = MattingSmooth(srcImg, smoothParam, cfg.shdwDegreeThres_matting, 
		cfg.maskErodeTimes_matting, cfg.maskDilateTimes_matting);
	if(cfg.debugImgsOutput)
	{
		ParamLoader::ShowParamInDir(cfg.debugImgsOutputDir, "step3_matting", smoothParam);
		cvi* resImg = RecoverProc::ImgRecoverNaive(srcImg, smoothParam, 1);
		cvsi(cfg.debugImgsOutputDir + "res_step3.png", resImg); cvri(resImg);
	}
	//cvri(smoothSrc);
	COutGreen("\r                                               \rStep 3: Matting res got.\n");

	return smoothParam;
}


cvi* DetectionProc::MRFSmooth(cvi* srcImg, cvi* paramResize, int nLabels)
{

	MRFProc proc;

	//cvCvtColor(srcImg, srcImg, CV_BGR2Lab);
	////MRF smoothing
	cvi* res;
	doFcvi(paramResize, i, j)
	{
		cvs2(paramResize, i, j, cvg2(paramResize, i, j));
	}
	proc.nChs = cfg.nCh_mrf; proc.colorSigma = cfg.colorsigma_mrfSTerm;
	proc.smoothWeight = cfg.sWeight_mrfSTerm;
	proc.SolveWithInitialAllCh(srcImg, paramResize, nLabels, res);
	cvCopy(res, paramResize);
	cvri(res);
	return paramResize;
	//cvCvtColor(srcImg, srcImg, CV_Lab2BGR);
	//cvri(paramResize);

// 	cvi* mask = cvci81(paramResize);
// 	doFcvi(mask, i, j)
// 	{
// 		cvs20(mask, i, j, cvg20(paramResize, i, j)*255);
// 	}
// 
// 	cvi* resMask;
// 	proc.SolveWithInitial(srcImg, NULL, mask, NULL, nLabels, resMask);
// 
// 	//cvSmooth(paramResize, paramResize, 1, srcImg->width / 30, srcImg->width / 30);
// 	doFcvi(mask, i, j)
// 	{
// 		auto v = cvg2(paramResize, i, j);
// 		v.val[0] = cvg20(resMask, i, j) / 255;
// 		cvs2(paramResize, i, j, v);
// 	}
// 	cvri(mask); cvri(resMask);
// 
// 	return paramResize;
}

cvi* DetectionProc::DecmpsSmooth(cvi* srcImg, cvi* param, int nLabels)
{
	MRFProc::guidInitMaxRatio_data = cfg.guidInitMaxRatio_data; 
	MRFProc::ratioSigma_data = cfg.ratioSigma_data;
	MRFProc::whiteWeight_data = cfg.whiteWeight_data;
	MRFProc::colorSigma_smooth = cfg.colorSigma_smooth;
	MRFProc::weight_smooth = cfg.weight_smooth;
// 	cout<<MRFProc::guidInitMaxRatio_data<<" "<<MRFProc::ratioSigma_data<<" "<<MRFProc::whiteWeight_data
// 		<<" "<<MRFProc::colorSigma_smooth<<" "<<MRFProc::weight_smooth;

	DecmpsProc proc;
	proc.cfg = cfg.dCfg;

	cvi* mask = cvci81(param);
	doFcvi(mask, i, j)
	{
		cvs20(mask, i, j, cvg20(param, i, j)*255);
	}

	cvi* resMask;
	proc.Analysis(srcImg, mask, resMask, cfg.dCfg.adjacentRatio, cfg.dCfg.nLabels);

	doFcvi(mask, i, j)
	{
		auto v = cvg2(param, i, j);
		v.val[0] = cvg20(resMask, i, j) / 255;
		if(v.val[0] > 0.95)
		{
			v.val[1] = 0;
			v.val[2] = 0;
		}
		cvs2(param, i, j, v);
	}
	cvri(mask); cvri(resMask);
	return param;

}

cvi* DetectionProc::MattingSmooth(cvi* srcImg, cvi* param, float shdwRatio, int erodeTimes, int dilateTimes)
{
	cvi* shdwMask = cvci81(param); cvZero(shdwMask);
	cvi* shdwMask2 = cvci(shdwMask);

	cvS avgShdw = cvs(0, 0, 0); int cnt = 0;

	//get shadow inner
	doFcvi(param, i, j)
	{
		if(cvg20(param, i, j) < shdwRatio)
		{
			cvs20(shdwMask, i, j, 255);
			avgShdw += cvg2(param, i, j); cnt++;
		}
		else if(cvg20(param, i, j) > 0.98) cvs20(shdwMask2, i, j, 255);
	}
	cvErode(shdwMask, shdwMask, 0, erodeTimes);
	cvErode(shdwMask2, shdwMask2, 0, dilateTimes);
	doFcvi(param, i, j)
	{
		if(cvg20(shdwMask2, i, j) == 255) cvs20(shdwMask, i, j, 0);
		else if(cvg20(shdwMask, i, j) == 0) cvs20(shdwMask, i, j, 128);
	}
	avgShdw /= cnt; //avgShdw = 0.45;

	cvsi("result//lp_src.png", srcImg);
	cvsi("result//lp_mask.png", shdwMask);

	//local propagation
	//image_editing::LocalEditpropProc proc;
	cvi* resMask = cvci(shdwMask);
	//proc.RunLocalEditprop(srcImg, shdwMask, resMask, 0.1f, 0.0001f, 0.1f);

	//alpha matting
	image_editing::AlphaMattingProc proc;
	cvi* back = cvci(srcImg), *fore = cvci(srcImg);
	proc.RunAlphaMatting(srcImg, shdwMask, back, fore, resMask);

	cvsi("result//lp_res.png", resMask);

	//update data
	doFcvi(resMask, i, j)
	{
		auto v = cvg2(param, i, j);
		auto v2 = cvg20(resMask, i, j) / 255;
		cvS v3 = avgShdw * v2 + cvs(1, 0, 0) * (1 - v2);
		//v.val[0] = v3.val[0];
		v = v3;
		if(v.val[0] > 0.999)
		{
			v.val[1] = 0;
			v.val[2] = 0;
		}
		cvs2(param, i, j, v);
	}
	cvri(shdwMask); cvri(resMask); cvri(shdwMask2);
	return param;
}

void DetectionProc::SetFileDir(string fileDir)
{
	//m_fileDir = fileDir;
	//m_srcDir = "01_VoteParam//"; // VoteParam PredictParam
	//m_resDir = "02_VoteDetection//";
	//wMkDir(m_fileDir + m_resDir);
}
