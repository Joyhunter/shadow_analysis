#include "StdAfx.h"
#include "gpm_analysis.h"

GPMAnalysisProc::GPMAnalysisProc(void)
{
}

GPMAnalysisProc::~GPMAnalysisProc(void)
{
}

void GPMAnalysisProc::SetFileDir(string fileDir)
{
	m_fileDir = fileDir;
}

void GPMAnalysisProc::ShdwAnlysis()
{
	wGetDirFiles(m_fileDir + "*.jpg", m_imgNames);

	param.range.setScale(0.67f, 1.5f);
	param.range.setRotate(-1.05f * (float)CV_PI, 1.05f * (float)CV_PI);
	param.range.setGain(cvs(0.2, 1, 1), cvs(5.0, 1, 1));
	param.range.setBias(cvs(0, -30, -30), cvs(0, 30, 30));
	param.downRatio = 6; param.multiLevelRatio = 2; param.extraLevelN = 0; 
	param.patchSize = 7;
	param.colorMode = CV_BGR2Lab;

	RunGPMForAllImages();
	
	GetStatistics();
}

string GPMAnalysisProc::GetCorrFileDir(string imgName)
{
	return "corr//" + imgName.substr(0, imgName.size()-4)
		+ "_" + toStr(param.downRatio) + "_" + toStr(param.multiLevelRatio) + "_" + toStr(param.extraLevelN)
		+ "_" + toStr(param.patchSize) + ".corr";
}
string GPMAnalysisProc::GetShadowSegDir(string imgName)
{
	return "seg//" + imgName.substr(0, imgName.size()-4) + "_shadow.png";
}
string GPMAnalysisProc::GetMaterialSegDir(string imgName)
{
	return "seg//" + imgName.substr(0, imgName.size()-4) + "_material.png";
}

void GPMAnalysisProc::RunGPMForAllImages()
{
	cout<<"Begin Running gpm for all images..."<<endl;
	doFv(i, m_imgNames)
	{
		cout<<"\r  Runing GPM for image "<<m_imgNames[i]<<"...\n";

		ImgContainer img(m_fileDir + m_imgNames[i], param.downRatio, param.colorMode);

		int gridSize = _i(img.src()->width * 0.28f);
		int gridOffset = _i(img.src()->width * 0.25f);
		
		string coorName = GetCorrFileDir(m_imgNames[i]);
		string corrFileName = m_fileDir + coorName;

		LmnIvrtPatchDistMetric metric2;
		GridGPMProc proc(&metric2, gridSize, gridOffset, 1, param.patchSize, 20, &(param.range));

		proc.RunGridGPMMultiScale(img, corrFileName, param.multiLevelRatio, param.extraLevelN);
		if(img.srcR() == NULL) img.GenerateResizedImg(1);
		proc.ShowGPMResUI(img, corrFileName, 80);

		cout<<"\r  GPM of "<<m_imgNames[i]<<" complete. File saved to "<<coorName<<".\n";
	}
}


void GPMAnalysisProc::GetStatistics()
{
	int distThres = 2000;

	GPMStatis sta;

	cout<<"Begin analysis matching result for all images...";
	doFv(i, m_imgNames)
	{
		GPMStatis sta1;

		string corrFile = GetCorrFileDir(m_imgNames[i]);
		string shadowMaskFile = GetShadowSegDir(m_imgNames[i]);
		string matMaskDir = GetMaterialSegDir(m_imgNames[i]);

		cvi* shMaskOri = cvlig(m_fileDir + shadowMaskFile);
		cvi* shMask = cvci81(shMaskOri->width / param.downRatio, shMaskOri->height / param.downRatio);
		cvResize(shMaskOri, shMask);
		cvi* matMaskOri = cvlig(m_fileDir + matMaskDir);
		cvi* matMask = cvci(shMask);
		cvResize(matMaskOri, matMask);
		cvri(shMaskOri); cvri(matMaskOri);

		DenseCorrBox2D box;
		box.Load(m_fileDir + corrFile);

		doFcvi(shMask, i, j)
		{
			MultiCorr corrs = box.GetCorrsPerGrid(i, j);

			int shdV = _i cvg20(shMask, i, j), matV = _i cvg20(matMask, i, j);
			doFv(k, corrs)
			{
				Corr& cor = corrs[k];
				if(cor.dist > distThres) continue;
				int shdV2 = _i cvg20(shMask, round(cor.x), round(cor.y));
				int matV2 = _i cvg20(matMask, round(cor.x), round(cor.y));
				if(shdV != shdV2 && matV == matV2) sta1.dssm++;
				else if(shdV != shdV2 && matV != matV2) sta1.dsdm++;
				else if(shdV == shdV2 && matV == matV2) sta1.sssm++;
				else sta1.ssdm++;
				sta1.matchN++;
			}
		}
		sta.add(sta1);

		cvri(shMask); cvri(matMask);
	}
	cout<<"\rAnalysis matching result for all images complete!"<<endl;
	sta.print();
}