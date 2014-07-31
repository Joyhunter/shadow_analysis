#include "StdAfx.h"
#include "shadow_analysis.h"

ShdwImgInfo Metric::AnalysisShadowFrame(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN)
{
	ShdwImgInfo info(histN);

	int patchRadius = 7;
	cvi* mask = cvci(imgShdwMask);
	cvErode(mask, mask, 0, patchRadius);

	int rsvN = sqr(2*patchRadius+1);
	vector<cvS> vShdw, vNonShdw;
	vShdw.reserve(rsvN); vNonShdw.reserve(rsvN);
	doFcvi(mask, i, j)
	{
		if(cvg20(mask, i, j) < 100) continue;
		info.nPixels++;

		//get pixels
		vShdw.resize(0); vNonShdw.resize(0);
		doFs(oi, -patchRadius, patchRadius) doFs(oj, -patchRadius, patchRadius)
		{
			if(!cvIn(i+oi, j+oj, mask)) continue;
			auto v1 = cvg2(imgShdw, i, j), v2 = cvg2(imgNonShdw, i, j);
			vShdw.push_back(v1); vNonShdw.push_back(v2);
		}

		ComputeInfoParam(vShdw, vNonShdw, info);
	}
	cvri(mask);
	return info;
}

//indoor 420828 pixels, 1.78778 assumption error
//outdoor 205145 pixels, 1.17416 assumption error
ShdwImgInfo RGBBiasMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN)
{
	return AnalysisShadowFrame(imgShdw, imgNonShdw, imgShdwMask, histN);
}

void RGBBiasMetric::ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info)
{
	//compute bias
	cvS sumShdw = cvs(0, 0, 0), sumUnShdw = cvs(0, 0, 0);
	doFv(i, vShdw)
	{
		sumShdw += vShdw[i]; sumUnShdw += vNonShdw[i];
	}
	bool swaped = false;
	if(cvSDSqr(sumShdw) > cvSDSqr(sumUnShdw))
	{
		swap(sumShdw, sumUnShdw);
		swaped = true;
	}
	vector<float> bias(3, 1.f);
	int histN = info.biasDstrbt[0].size();
	doF(k, 3)
	{
		if(sumUnShdw.val[k] != 0)
		{
			bias[k] = _f sumShdw.val[k] / _f sumUnShdw.val[k];
		}
		int idx = clamp(_i floor(bias[k]*histN), 0, histN-1);
		info.biasDstrbt[k][idx]++;
		if(swaped) bias[k] = (bias[k] == 0)?255:(1.f/bias[k]); 
	}

	//get assumption error
	doF(k, 3)
	{
		if(bias[k] < 1e-3f) bias[k] = 1e-3f;
	}
	double dist = 0;
	cvS bia = cvs(bias[0], bias[1], bias[2]);
	doFv(k, vShdw)
	{
		dist += cvSD(clamp(vShdw[k] / bia, cvs(0, 0, 0), cvs(255, 255, 255)), vNonShdw[k]);
	}
	dist = (vShdw.size() == 0)?0:dist/vShdw.size();
	info.asmptDist += dist;
}

ShdwImgInfo HLSGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN)
{
	cvi* i1 = cvci(imgShdw), *i2 = cvci(imgNonShdw);
	cvCvtColor(imgShdw, i1, CV_BGR2HLS_FULL);
	cvCvtColor(imgNonShdw, i2, CV_BGR2HLS_FULL);
	ShdwImgInfo res = AnalysisShadowFrame(i1, i2, imgShdwMask, histN);
	cvri(i1); cvri(i2);
	return res;
}

void HLSGainMetric::ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info)
{
	//compute gain
	cvS sumShdw = cvs(0, 0, 0), sumUnShdw = cvs(0, 0, 0);
	doFv(i, vShdw)
	{
		sumShdw += vShdw[i]; sumUnShdw += vNonShdw[i];
	}
	bool swaped = false;
	if(cvSDSqr(sumShdw) > cvSDSqr(sumUnShdw))
	{
		swap(sumShdw, sumUnShdw);
		swaped = true;
	}
	vector<float> bias(3, 1.f);
	int histN = info.biasDstrbt[0].size();
	doF(k, 3)
	{
		if(sumUnShdw.val[k] != 0)
		{
			bias[k] = _f sumShdw.val[k] / _f sumUnShdw.val[k];
		}
		int idx = clamp(_i floor(bias[k]*histN), 0, histN-1);
		info.biasDstrbt[k][idx]++;
		if(swaped) bias[k] = (bias[k] == 0)?255:(1.f/bias[k]); 
	}

	//get assumption error
	doF(k, 3)
	{
		if(bias[k] < 1e-3f) bias[k] = 1e-3f;
	}
	double dist = 0;
	cvS bia = cvs(bias[0], bias[1], bias[2]);
	doFv(k, vShdw)
	{
		dist += cvSD(clamp(vShdw[k] / bia, cvs(0, 0, 0), cvs(255, 255, 255)), vNonShdw[k]);
	}
	dist = (vShdw.size() == 0)?0:dist/vShdw.size();
	info.asmptDist += dist;
}

ShdwImgInfo RGBBiasGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN)
{
	return AnalysisShadowFrame(imgShdw, imgNonShdw, imgShdwMask, histN);
}

void RGBBiasGainMetric::ComputeInfoParam(vector<cvS>& vShdw, vector<cvS>& vNonShdw, ShdwImgInfo& info)
{

}

ShdwAnlysisProc::ShdwAnlysisProc(void)
{
}

ShdwAnlysisProc::~ShdwAnlysisProc(void)
{
}

void ShdwAnlysisProc::SetFileDir(string fileDir)
{
	m_fileDir = fileDir;
}

void ShdwAnlysisProc::ShdwAnlysis()
{
	wGetDirFiles(m_fileDir + "*_n.jpg", m_imgNames);
	
	RGBBiasMetric metric1;
	HLSGainMetric metric2;
	RGBBiasGainMetric metric3;
	int histN = 128;

	ShdwImgInfo sumInfo(histN);
	vector<ShdwImgInfo> info(m_imgNames.size(), histN);

	omp_set_num_threads(1);
#pragma omp parallel for
	doF(i, _i m_imgNames.size())
	{
		cout<<"\rProcessing "<<m_imgNames[i]<<"..";
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		string imgShdwDir = m_fileDir + imgPrefix + ".jpg";
		string imgUnShdwDir = m_fileDir + imgPrefix + "_n.jpg";
		string imgShdwMaskDir = m_fileDir + "mask//" + imgPrefix + "_mask.png";

		cvi* imgShdw = cvlic(imgShdwDir);
		cvi* imgNonShdw = cvlic(imgUnShdwDir);
		cvi* imgShdwMask = cvlig(imgShdwMaskDir);
		float resizeRatio = 0.4f;
		cvi* i1 = cvci83(_i (imgShdw->width * resizeRatio), _i (imgShdw->height * resizeRatio));
		cvi* i2 = cvci83(i1), *i3 = cvci81(i1);
		cvResize(imgShdw, i1); cvResize(imgNonShdw, i2); cvResize(imgShdwMask, i3);

		info[i] = metric1.AnalysisShadow(i1, i2, i3, histN);
		
		cvi* res = VislzDstrbt(info[i].biasDstrbt, info[i].nPixels);
		cvsi("result//" + imgPrefix + ".png", res);
		cvri(res);

		cvri(imgShdw); cvri(imgNonShdw); cvri(imgShdwMask); cvri(i1); cvri(i2); cvri(i3);
		//cout<<endl<<m_imgNames[i]<<" "<<info[i].nPixels<<" "<<info[i].asmptDist<<endl;
	}
	doFv(i, info) sumInfo += info[i];

	cvi* res = VislzDstrbt(sumInfo.biasDstrbt, sumInfo.nPixels);
	cvsi("result//sum.png", res);
	cvri(res);
	cout<<sumInfo.nPixels<<" "<<sumInfo.asmptDist<<endl;

	cout<<"\rDone.";
}

cvi* ShdwAnlysisProc::VislzDstrbt(vector<vector<int> >& dstrbt, int nPixels)
{
	int chN = dstrbt.size();
	int histN = dstrbt[0].size();

	int histWidth = max2(200 / histN, 3);
	int histHeight = 100;
	int showHRatio = histN / 20;
	int gap = 1;
	
	cvi* img = cvci83(histWidth * histN, histHeight * chN);
	cvZero(img);

	doF(k, chN)
	{
		doF(i, histN)
		{
			cvS color = cvs(0, 0, 255);
			if(k == 1) color = cvs(0, 255, 0);
			if(k == 2) color = cvs(255, 0, 0);
			color = color * i / histN;
			int h = clamp(showHRatio * dstrbt[k][i] * histHeight / nPixels, 0, histHeight-gap);
			cvRectangle(img, cvPoint(i*histWidth + gap, histHeight-h+k*histHeight),
				cvPoint((i+1)*histWidth - gap, histHeight+k*histHeight), color, -1);
		}
	}
	return img;
}