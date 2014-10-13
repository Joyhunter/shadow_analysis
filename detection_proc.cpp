#include "StdAfx.h"
#include "detection_proc.h"
#include "param_loader.h"
#include "../shadow_removal/mrf_proc.h"
#include "../shadow_removal/decmps_proc.h"

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
	m_resDir = "02_Detection//";
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
		if(imgPrefix != "004") continue;
		cout<<"Handling "<<imgPrefix<<".\n";

		//load source image and param prediction
		cvi* srcImg = cvlic(m_fileDir + imgPrefix + ".png");
		cvi* param = ParamLoader::LoadParamFromDir(m_fileDir + m_srcDir, imgPrefix);

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
	cvi* smoothParam = MRFSmooth(srcImg, param, 16);
	smoothParam = DecmpsSmooth(srcImg, smoothParam, 128);

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

	cvSmooth(paramResize, paramResize, 1, srcImg->width / 30, srcImg->width / 30);
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