#include "StdAfx.h"
#include "shadow_analysis.h"

ShdwImgInfo Metric::AnalysisShadowBase(const cvi* _imgShdw, const cvi* _imgNonShdw, const cvi* imgShdwMask, int histN, 
	int patchRadius, vector<bool>& useGain, vector<bool>& useBias, int cSpaceType, int cSpaceType2)
{
	cvi* imgShdw = cvci(_imgShdw), *imgNonShdw = cvci(_imgNonShdw);
	cvi* imgShdwRes = cvci(_imgShdw), *imgNonShdwRes = cvci(_imgNonShdw);
	if(cSpaceType > 0)
	{
		cvCvtColor(_imgShdw, imgShdw, cSpaceType);
		cvCvtColor(_imgNonShdw, imgNonShdw, cSpaceType);
	}

	ShdwImgInfo info(histN);

	cvi* mask = cvci(imgShdwMask);
	cvErode(mask, mask, 0, patchRadius);

	int rsvN = sqr(2*patchRadius+1);

	doFcvi(mask, i, j)
	{
		if(cvg20(mask, i, j) < 100) continue;
		info.nPixels++;

		//get pixels
		vector<cvS> vShdw, vNonShdw;
		vShdw.reserve(rsvN); vNonShdw.reserve(rsvN);
		int curIdx = 0;
		doFs(oi, -patchRadius, patchRadius) doFs(oj, -patchRadius, patchRadius)
		{
			if(!cvIn(i+oi, j+oj, mask)) continue;
			if(oi == 0 && oj == 0) curIdx = vShdw.size();
			auto v1 = cvg2(imgShdw, i+oi, j+oj), v2 = cvg2(imgNonShdw, i+oi, j+oj);
			vShdw.push_back(v1); vNonShdw.push_back(v2);
		}

		ComputeInfoParam(vShdw, vNonShdw, info, useGain, useBias, cSpaceType2);

		cvs2(imgShdwRes, i, j, vShdw[curIdx]);
		cvs2(imgNonShdwRes, i, j, vNonShdw[curIdx]);
	}

// 	cvsi("_shdw_before.png", _imgShdw);
// 	cvsi("_recover1.png", imgShdwRes);
// 	cvsi("_nonshdw_before.png", _imgNonShdw);
// 	cvsi("_recover2.png", imgNonShdwRes);
// 	pause;

	cvri(mask);
	cvri(imgShdw); cvri(imgNonShdw);
	cvri(imgShdwRes); cvri(imgNonShdwRes);
	info.asmptDist /= info.nPixels;
	return info;
}

//bias && gain : unshadow = (shadow - avgSh) / gain + avgSh + bias;  gain: 0.1~1~10 bias: 0~255
//!bias&& gain : unshadow = shadow / gain; gain: 0~1
//bias &&!gain : unshadow = shadow + bias; bias: -255~255
//!bias&&!gain : unshadow = shadow;
void Metric::ComputeInfoParam(vector<cvS>& _vShdw, vector<cvS>& _vNonShdw, ShdwImgInfo& info,
	vector<bool>& useGain, vector<bool>& useBias, int cSpaceType2)
{
	int N = _vShdw.size();
	vector<cvS>* vShdw = &_vShdw, *vNonShdw = &_vNonShdw;

	//compute avg, sigma
	cvS avgShdw = cvs(0, 0, 0), avgUnShdw = cvs(0, 0, 0);
	cvS sigShdw = cvs(0, 0, 0), sigUnShdw = cvs(0, 0, 0);
	doFv(i, (*vShdw))
	{
		avgShdw += (*vShdw)[i]; avgUnShdw += (*vNonShdw)[i];
	}
	avgShdw /= N; avgUnShdw /= N;
// 	if(cvSD(avgShdw) > cvSD(avgUnShdw))
// 	{
// 		vShdw = &_vNonShdw;
// 		vNonShdw = &_vShdw;
// 		swap(avgShdw, avgUnShdw);
// 	}
	doFv(i, (*vShdw))
	{
		sigShdw += sqr((*vShdw)[i] - avgShdw);
		sigUnShdw += sqr((*vNonShdw)[i] - avgUnShdw);
	}
	sigShdw = (sigShdw / N) ^ 0.5;
	sigUnShdw = (sigUnShdw / N) ^ 0.5;

	//get bias, gain    
	double maxGain = 10;
	double maxBias = 255;
	int histN = info.biasDstrbt[0].size();
	vector<double> biasV(3), gainV(3);
	doF(k, 3)
	{
		double gain = 1, bias = 0;
		int gainIdx = -1, biasIdx = -1;
		if(useGain[k] && useBias[k])
		{
			if(sigUnShdw.val[k] > 1e-6) gain = sigShdw.val[k] / sigUnShdw.val[k];
			gain = clamp(gain, 1.0 / maxGain, maxGain);
			bias = avgUnShdw.val[k] - avgShdw.val[k]; //cout<<gain<<endl; pause;
			gainIdx = clamp(_i floor((log(gain) + log(maxGain)) / 2 / log(maxGain) * histN), 0, histN-1);
			biasIdx = clamp(_i floor(bias*histN/maxBias), 0, histN-1);
		}
		if(useGain[k] && !useBias[k])
		{
			gain = avgShdw.val[k] / avgUnShdw.val[k];
			gainIdx = clamp(_i floor(gain*histN), 0, histN-1);
			//gainIdx = clamp(_i floor((log(gain) + log(maxGain)) / 2 / log(maxGain) * histN), 0, histN-1);
		}
		if(!useGain[k] && useBias[k])
		{
			bias = avgUnShdw.val[k] - avgShdw.val[k];
			biasIdx = clamp(_i floor((bias+maxBias)*histN/maxBias/2), 0, histN-1);
		}
		if(biasIdx >= 0) info.biasDstrbt[k][biasIdx]++;
		if(gainIdx >= 0) info.gainDstrbt[k][gainIdx]++;
		biasV[k] = bias; gainV[k] = gain;
	}

	//get assumption error
	vector<cvS> vRecUnShadow((*vShdw).size());
	doFv(i, (*vShdw))
	{
		vector<double> temp(3);
		doF(k, 3)
		{
			if(useGain[k] && useBias[k])
			{
				temp[k] = clamp(((*vShdw)[i].val[k]-avgShdw.val[k]) / max2(gainV[k], 1e-3) 
					+ avgShdw.val[k] + biasV[k], 0.0, 255.0);
			}
			if(useGain[k] && !useBias[k])
			{
				temp[k] = clamp((*vShdw)[i].val[k] / gainV[k], 0.0, 255.0);
			}
			if(!useGain[k] && useBias[k])
			{
				temp[k] = clamp((*vShdw)[i].val[k] + biasV[k], 0.0, 255.0);
			}
			if(!useGain[k] && !useBias[k])
			{
				temp[k] = (*vShdw)[i].val[k];
			}
		}
		vRecUnShadow[i] = cvs(temp[0], temp[1], temp[2]);
	}

	if(cSpaceType2 > 0)
	{
		cvi* temp = cvci83(2, (*vShdw).size());
		doFv(i, (*vShdw))
		{
			cvs2(temp, i, 0, vRecUnShadow[i]);
			cvs2(temp, i, 1, (*vNonShdw)[i]);
		}
		cvCvtColor(temp, temp, cSpaceType2);
		doFv(i, (*vShdw))
		{
			vRecUnShadow[i] = cvg2(temp, i, 0);
			(*vNonShdw)[i] = cvg2(temp, i, 1);
		}
	}

	double dist = 0;
	doFv(i, (*vShdw))
	{
		dist += cvSD(vRecUnShadow[i], (*vNonShdw)[i]);
		(*vShdw)[i] = vRecUnShadow[i];
	}
	dist = ((*vShdw).size() == 0)?0:dist/(*vShdw).size();
	info.asmptDist += dist;
}

//indoor 420828 pixels, 8.1612 assumption error
//outdoor 205145 pixels, 12.8593 assumption error
//R(gain)G(gain)B(gain)
ShdwImgInfo RGBGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN, int patchRadius)
{
	vector<bool> useGain(3, true), useBias(3, false); 
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, -1, -1);
}

//indoor 9.37166
//outdoor 14.9699 15.9736(without L bias)
//H(gain bias)L(bias)S(bias)
ShdwImgInfo HLSGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN, int patchRadius)
{
	vector<bool> useGain(3, false), useBias(3, true);
	useGain[1] = true; //useBias[0] = false; useBias[2] = false;
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, 
		CV_BGR2HLS_FULL, CV_HLS2BGR_FULL);
}

//indoor 420828 pixels, 7.43204 assumption error
//outdoor 205145 pixels, 12.0348 assumption error
//R(gain, bias)G(gain, bias)B(gain, bias)
ShdwImgInfo RGBBiasGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN, int patchRadius)
{
	vector<bool> useGain(3, true), useBias(3, true); 
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, -1, -1);
}

//indoor 8.52214(without L bias)
//outdoor 12.23 13.3504(without L bias)
//L(gain bias)A(bias)B(bias)
ShdwImgInfo LABGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, int histN, int patchRadius)
{
	vector<bool> useGain(3, false), useBias(3, false);
	useGain[0] = true; useBias[1] = true; useBias[2] = true;
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, 
		CV_BGR2Lab, CV_Lab2BGR);
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
	
	RGBGainMetric metric1;
	HLSGainMetric metric2;
	RGBBiasGainMetric metric3;
	LABGainMetric metric4;
	Metric* metric = &metric4;
	int histN = 128;
	float resizeRatio = 0.4f;
	int patchRadius = 7;

	ShdwImgInfo sumInfo(histN);
	vector<ShdwImgInfo> info(m_imgNames.size(), histN);

	omp_set_num_threads(1);
#pragma omp parallel for
	doF(i, _i m_imgNames.size())
	{
		//cout<<"\rProcessing "<<m_imgNames[i]<<"..";
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		string imgShdwDir = m_fileDir + imgPrefix + ".jpg";
		string imgUnShdwDir = m_fileDir + imgPrefix + "_n.jpg";
		string imgShdwMaskDir = m_fileDir + "mask//" + imgPrefix + "_mask.png";

		cvi* imgShdw = cvlic(imgShdwDir);
		cvi* imgNonShdw = cvlic(imgUnShdwDir);
		cvi* imgShdwMask = cvlig(imgShdwMaskDir);
		cvi* i1 = cvci83(_i (imgShdw->width * resizeRatio), _i (imgShdw->height * resizeRatio));
		cvi* i2 = cvci83(i1), *i3 = cvci81(i1);
		cvResize(imgShdw, i1); cvResize(imgNonShdw, i2); cvResize(imgShdwMask, i3);

		info[i] = metric->AnalysisShadow(i1, i2, i3, histN, patchRadius);
		//info[i] = metric3.AnalysisShadow(imgShdw, imgNonShdw, imgShdwMask, histN);

		cvi* res = VislzDstrbt(info[i].biasDstrbt, info[i].nPixels);
		cvsi("result//" + imgPrefix + "_bias.png", res);
		cvri(res);
		cvi* res2 = VislzDstrbt(info[i].gainDstrbt, info[i].nPixels);
		cvsi("result//" + imgPrefix + "_gain.png", res2);
		cvri(res2);

		cvri(imgShdw); cvri(imgNonShdw); cvri(imgShdwMask); cvri(i1); cvri(i2); cvri(i3);
		cout<<m_imgNames[i]<<": "<<info[i].nPixels<<" pixels, avg Patch error = "<<info[i].asmptDist<<"."<<endl;
	}
	doFv(i, info) sumInfo += info[i];

	cvi* res = VislzDstrbt(sumInfo.biasDstrbt, sumInfo.nPixels);
	cvsi("result//sum_bias.png", res);
	cvri(res);
	cvi* res2 = VislzDstrbt(sumInfo.gainDstrbt, sumInfo.nPixels);
	cvsi("result//sum_gain.png", res2);
	cvri(res2);
	cout<<"Avg patch error per image: "<<sumInfo.asmptDist/info.size()<<endl;

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
			//color = color * i / histN;
			int h = clamp(showHRatio * dstrbt[k][i] * histHeight / nPixels, 0, histHeight-gap);
			cvRectangle(img, cvPoint(i*histWidth + gap, histHeight-h+k*histHeight),
				cvPoint((i+1)*histWidth - gap, histHeight+k*histHeight), color, -1);
		}
	}
	return img;
}