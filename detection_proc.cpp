#include "StdAfx.h"
#include "detection_proc.h"
#include "param_loader.h"
#include "../shadow_removal/mrf_proc.h"
#include "../shadow_removal/decmps_proc.h"
#include "recovery_proc.h"
#include "local_propagation/local_editprop_proc.h"
#include "alpha_matting/alpha_matting_proc.h"

DetectionProc::DetectionProc(void)
{
}


DetectionProc::~DetectionProc(void)
{
}

void DetectionProc::SetFileDir(string fileDir)
{
	m_fileDir = fileDir;
	m_srcDir = "01_PredictParam//"; // VoteParam PredictParam
	m_resDir = "02_VoteDetection//";
	wMkDir(m_fileDir + m_resDir);
}

void DetectionProc::DetectCastShadow()
{
	wGetDirFiles(m_fileDir + "*_n.png", m_imgNames);

	doFv(i, m_imgNames)
	{
		//get image prefix
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		if(imgPrefix != "000") continue;
		cout<<"Handling "<<imgPrefix<<".\n";

		//load source image and param prediction
		cvi* srcImg = cvlic(m_fileDir + imgPrefix + ".png");
		cvi* param = ParamLoader::LoadParamFromDir(m_fileDir + m_srcDir, imgPrefix);//m_srcDir

		//cast shadow detect
		cvi* paramNew = DetectCastShadow(srcImg, param);

		//save detect result
		ParamLoader::SaveParamToDir(m_fileDir + m_resDir, imgPrefix, paramNew);
		ParamLoader::ShowParamInDir(m_fileDir + m_resDir, imgPrefix, paramNew);

		cvri(srcImg); cvri(param); cvri(paramNew);
	}
}

cvi* DetectionProc::DetectCastShadow(cvi* srcImg, cvi* param)
{

	//srcImg = cvlic("result//004blur2.png");

	//cvi* smoothParam = cvci(param);
	cvi* smoothParam = MRFSmooth(srcImg, param, 16);

	//Test
	cvi* resImg = RecoverProc::ImgRecoverNaive(srcImg, smoothParam, 1);
	cvsi("result//1.png", resImg);
	
	//smoothing
	smoothParam = DecmpsSmooth(srcImg, smoothParam, 16);

	//Test
	resImg = RecoverProc::ImgRecoverNaive(srcImg, smoothParam, 1);
	cvsi("result//2.png", resImg); cvri(resImg);

	//cvi* smoothSrc = cvci(srcImg);
	//cvSmooth(srcImg, smoothSrc, 2, 31, 31);
	smoothParam = MattingSmooth(srcImg, smoothParam, 0.8f, 5, 10);

	//Test
	resImg = RecoverProc::ImgRecoverNaive(srcImg, smoothParam, 1);
	cvsi("result//3.png", resImg); cvri(resImg);

	//cvri(smoothSrc);

	return smoothParam;
}


cvi* DetectionProc::MRFSmooth(cvi* srcImg, cvi* param, int nLabels)
{
	cvi* paramResize = cvci323(srcImg);
	cvResize(param, paramResize);

	MRFProc proc;

	//cvCvtColor(srcImg, srcImg, CV_BGR2Lab);
	//cvi* res;
	////MRF smoothing
	//proc.SolveWithInitialAllCh(srcImg, paramResize, 64, res);
	//cvCvtColor(srcImg, srcImg, CV_Lab2BGR);
	//cvri(paramResize);

	cvi* mask = cvci81(paramResize);
	doFcvi(mask, i, j)
	{
		cvs20(mask, i, j, cvg20(paramResize, i, j)*255);
	}

	cvi* resMask;
	proc.SolveWithInitial(srcImg, NULL, mask, NULL, nLabels, resMask);

	//cvSmooth(paramResize, paramResize, 1, srcImg->width / 30, srcImg->width / 30);
	doFcvi(mask, i, j)
	{
		auto v = cvg2(paramResize, i, j);
		v.val[0] = cvg20(resMask, i, j) / 255;
		cvs2(paramResize, i, j, v);
	}
	cvri(mask); cvri(resMask);

	return paramResize;
}

cvi* DetectionProc::DecmpsSmooth(cvi* srcImg, cvi* param, int nLabels)
{
	DecmpsProc proc;

	cvi* mask = cvci81(param);
	doFcvi(mask, i, j)
	{
		cvs20(mask, i, j, cvg20(param, i, j)*255);
	}

	cvi* resMask;
	proc.Analysis(srcImg, mask, resMask, 1.0f, nLabels);

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