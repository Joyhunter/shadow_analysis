#include "StdAfx.h"
#include "shadow_analysis.h"
#include "param_loader.h"

//---------------------- ShdwImgInfo ----------------------------

//bias && gain : unshadow = (shadow - avgSh) / gain + avgSh + bias;  gain: 0.1~1~10 bias: 0~255
//!bias&& gain : unshadow = shadow / gain; gain: 0~1
//bias &&!gain : unshadow = shadow + bias; bias: -255~255
//!bias&&!gain : unshadow = shadow;
double maxBias = 255;
double maxGain = 10;
int ShdwImgInfo::computeBiasIdx(double bias, int histN, bool useBias, bool useGain)
{
	//0~255
	if(useBias && useGain)
		return clamp(_i floor(bias*histN/maxBias), 0, histN-1);
	//-255~255
	else if(useBias && !useGain)
		return clamp(_i floor((bias+maxBias)*histN/maxBias/2), 0, histN-1);
	return -1;
}
int ShdwImgInfo::computeGainIdx(double gain, int histN, bool useBias, bool useGain)
{
	//0.1~1~10
	if(useBias && useGain)
		return clamp(_i floor((log(gain) + log(maxGain)) / 2 / log(maxGain) * histN), 0, histN-1);
	//0~1
	else if(!useBias && useGain)
		return clamp(_i floor(gain*histN), 0, histN-1);
	return -1;
}
double ShdwImgInfo::computeBiasV(int biasIdx, int histN, bool useBias, bool useGain)
{
	double maxBias = 255;
	if(useBias && useGain)
		return maxBias / histN * (_d biasIdx + 0.5);
	else if(useBias && !useGain)
		return maxBias * 2 / histN * (_d biasIdx + 0.5) - maxBias;
	return 0;
}
double ShdwImgInfo::computeGainV(int gainIdx, int histN, bool useBias, bool useGain)
{
	if(useBias && useGain)
		return exp((_d gainIdx + 0.5) / histN * 2 * log(maxGain) - log(maxGain));
	else if(!useBias && useGain)
		return (_d gainIdx + 0.5) / histN;
	return 1;
}

void getMinMaxIdx(vector<int> distri, int N, float minMaxRatio, int& minIdx, int& maxIdx)
{
	minIdx = -1;
	int val = 0;
	while(minIdx < _i distri.size()-1 && val < _f N * minMaxRatio)
	{
		minIdx++;
		val += distri[minIdx];
	}
	maxIdx = distri.size(); val = 0;
	while(maxIdx >= 1 && val < _f N * minMaxRatio)
	{
		maxIdx--;
		val += distri[maxIdx];
	}
}
void ShdwImgInfo::computeMinMax(float minMaxRatio)
{
	biasMin.resize(3, -1); biasMax.resize(3, -1);
	gainMin.resize(3, -1); gainMax.resize(3, -1);
	doF(k, 3)
	{
		if(useBias[k])
		{
			int minIdx, maxIdx;
			getMinMaxIdx(biasDstrbt[k], nPixels, minMaxRatio, minIdx, maxIdx); //cout<<minIdx<<" "<<maxIdx<<endl;
			biasMin[k] = computeBiasV(minIdx, biasDstrbt[k].size(), true, useGain[k]);
			biasMax[k] = computeBiasV(maxIdx, biasDstrbt[k].size(), true, useGain[k]);
		}
		if(useGain[k])
		{
			int minIdx, maxIdx;
			getMinMaxIdx(gainDstrbt[k], nPixels, minMaxRatio, minIdx, maxIdx); //cout<<minIdx<<" "<<maxIdx<<endl;
			gainMin[k] = computeGainV(minIdx, gainDstrbt[k].size(), useBias[k], true);
			gainMax[k] = computeGainV(maxIdx, gainDstrbt[k].size(), useBias[k], true);
		}
	}
}
void ShdwImgInfo::outputMinMax()
{
	doF(k, 3)
	{
		if(useBias[k]) cout<<"Channel "<<k<<" bias: ("<<biasMin[k]<<", "<<biasMax[k]<<").\n";
		if(useGain[k]) cout<<"Channel "<<k<<" gain: ("<<gainMin[k]<<", "<<gainMax[k]<<").\n";
	}
}

//---------------------- RGBRatioConverter ----------------------------

void RGBRatioConverter::convert(const cvi* src, cvi* dst)
{
	doFcvi(src, i, j)
	{
		cvS v = cvg2(src, i, j);
		cvS v2 = v;
		v2.val[0] = (v.val[0] + v.val[1] + v.val[2]) / 3;
		double base = max2(v.val[0], 1.0);
		v2.val[1] = (clamp(log(v.val[1] / base) / log(10.0), -1.0, 1.0) * 127.5 + 127.5) / 2;
		v2.val[2] = (clamp(log(v.val[2] / base) / log(10.0), -1.0, 1.0) * 127.5 + 127.5) / 2;
		cvs2(dst, i, j, v2);
	}
}
void RGBRatioConverter::invConvert(const cvi* src, cvi* dst)
{
	doFcvi(src, i, j)
	{
		cvS v = cvg2(src, i, j);
		cvS v2 = v;

		double sumA = v.val[0] * 3;
		double ratio1 = exp((v.val[1] - 127.5) / 127.5 * log(10.0));
		double ratio2 = exp((v.val[2] - 127.5) / 127.5 * log(10.0));

		v2.val[0] = sumA / (1+ratio1+ratio2);
		v2.val[1] = v2.val[0] * ratio1;
		v2.val[2] = v2.val[0] * ratio2;
		v2 = clamp(v2, cvs(0, 0, 0), cvs(255, 255, 255));

		cvs2(dst, i, j, v2);
	}
}

//---------------------- Metric ---------------------------------------

ShdwImgInfo Metric::AnalysisShadowBase(const cvi* _imgShdw, const cvi* _imgNonShdw, const cvi* imgShdwMask, int histN, 
	int patchRadius, vector<bool>& useGain, vector<bool>& useBias, int cSpaceType, int cSpaceType2, 
	OUT vector<double>& gbVs, ImgConverter* converter)
{
	cvi* imgShdw = cvci(_imgShdw), *imgNonShdw = cvci(_imgNonShdw);
	cvi* imgShdwRes = cvci(_imgShdw), *imgNonShdwRes = cvci(_imgNonShdw);
	if(cSpaceType > 0)
	{
		cvCvtColor(_imgShdw, imgShdw, cSpaceType);
		cvCvtColor(_imgNonShdw, imgNonShdw, cSpaceType);
	}
	else if(converter != NULL)
	{
		converter->convert(_imgShdw, imgShdw);
		converter->convert(_imgNonShdw, imgNonShdw);
	}

	ShdwImgInfo info(histN);

	cvi* mask = cvci(imgShdwMask);
	cvDilate(mask, mask, 0, patchRadius);
	//cvErode(mask, mask, 0, patchRadius);

	int rsvN = sqr(2*patchRadius+1);

	cvi* bgShowImg[6];
	doF(k, 6)
	{
		bgShowImg[k] = cvci(mask);
		cvZero(bgShowImg[k]);
	}
	doFcvi(mask, i, j)
	{
		if(cvg20(mask, i, j) > 128)
		{
// 			doF(k, 3)
// 			{
// 				gbVs[(i*mask->width+j)*6+2*k] = 1;
// 				gbVs[(i*mask->width+j)*6+2*k+1] = 0;
// 			}
// 			continue;
		}
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

		vector<int> idxs(6, -1);
		vector<double> gbV(6, 0);
		ComputeInfoParam(vShdw, vNonShdw, info, useGain, useBias, cSpaceType2, idxs, gbV, converter);
		doF(k, 6)
		{
			gbVs[(i*mask->width+j)*6+k] = gbV[k]; 
		}
		
		doF(k, 6)
		{
			if(idxs[k] != -1) cvs20(bgShowImg[k], i, j, idxs[k]*255/histN);
		}

		cvs2(imgShdwRes, i, j, vShdw[curIdx]);
		cvs2(imgNonShdwRes, i, j, vNonShdw[curIdx]);
	}

// 	cvsi("biasGain_src.png", imgShdw);
// 	doF(k, 6)
// 	{
// 		cvsi("biasGain_" + toStr(k) + ".png", bgShowImg[k]);
// 		cvri(bgShowImg[k]);
// 	}pause;

// 	cvsi("_shdw_before.png", _imgShdw);
// 	cvsi("_recover1.png", imgShdwRes);
// 	cvsi("_nonshdw_before.png", _imgNonShdw);
// 	cvsi("_recover2.png", imgNonShdwRes);
// 	pause;

	cvri(mask);
	cvri(imgShdw); cvri(imgNonShdw);
	cvri(imgShdwRes); cvri(imgNonShdwRes);
	info.asmptDist /= info.nPixels;
	info.useBias = useBias;
	info.useGain = useGain;
	return info;
}

void Metric::ComputeInfoParam(vector<cvS>& _vShdw, vector<cvS>& _vNonShdw, ShdwImgInfo& info,
	vector<bool>& useGain, vector<bool>& useBias, int cSpaceType2, vector<int>& idxs, vector<double>& gbV, 
	ImgConverter* converter)
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
			gainIdx = info.computeGainIdx(gain, histN, true, true);
			//clamp(_i floor((log(gain) + log(maxGain)) / 2 / log(maxGain) * histN), 0, histN-1);
			biasIdx = info.computeBiasIdx(bias, histN, true, true);
			//clamp(_i floor(bias*histN/maxBias), 0, histN-1);
		}
		if(useGain[k] && !useBias[k])
		{
			gain = avgShdw.val[k] / avgUnShdw.val[k];
			gainIdx = info.computeGainIdx(gain, histN, false, true);
			//clamp(_i floor(gain*histN), 0, histN-1);
			//gainIdx = clamp(_i floor((log(gain) + log(maxGain)) / 2 / log(maxGain) * histN), 0, histN-1);
		}
		if(!useGain[k] && useBias[k])
		{
			bias = avgUnShdw.val[k] - avgShdw.val[k];
			biasIdx = info.computeBiasIdx(bias, histN, true, false);
			//clamp(_i floor((bias+maxBias)*histN/maxBias/2), 0, histN-1);
		}
		if(biasIdx >= 0) info.biasDstrbt[k][biasIdx]++;
		if(gainIdx >= 0) info.gainDstrbt[k][gainIdx]++;
		biasV[k] = bias; gainV[k] = gain;
		idxs[2*k] = biasIdx; idxs[2*k+1] = gainIdx; 
		gbV[2*k] = gain; gbV[2*k+1] = bias;
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

	if(cSpaceType2 > 0 || converter != NULL)
	{
		cvi* temp = cvci83(2, (*vShdw).size());
		doFv(i, (*vShdw))
		{
			cvs2(temp, i, 0, vRecUnShadow[i]);
			cvs2(temp, i, 1, (*vNonShdw)[i]);
		}
		if(cSpaceType2 > 0) cvCvtColor(temp, temp, cSpaceType2);
		else converter->invConvert(temp, temp);
		doFv(i, (*vShdw))
		{
			vRecUnShadow[i] = cvg2(temp, i, 0);
			(*vNonShdw)[i] = cvg2(temp, i, 1);
		}
		cvri(temp);
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

ShdwImgInfo RGBGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
	OUT vector<double>& gbVs, int histN, int patchRadius)
{
	vector<bool> useGain(3, true), useBias(3, false); 
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, -1, -1, gbVs);
}

ShdwImgInfo HLSGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
	OUT vector<double>& gbVs, int histN, int patchRadius)
{
	vector<bool> useGain(3, false), useBias(3, true);
	useGain[1] = true; //useBias[0] = false; useBias[2] = false;
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, 
		CV_BGR2HLS_FULL, CV_HLS2BGR_FULL, gbVs);
}

ShdwImgInfo RGBBiasGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
	OUT vector<double>& gbVs, int histN, int patchRadius)
{
	vector<bool> useGain(3, true), useBias(3, true); 
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, -1, -1, gbVs);
}

ShdwImgInfo LABGainMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
	OUT vector<double>& gbVs, int histN, int patchRadius)
{
	vector<bool> useGain(3, false), useBias(3, false);
	useGain[0] = true; useBias[1] = true; useBias[2] = true;
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, 
		CV_BGR2Lab, CV_Lab2BGR, gbVs);
}

ShdwImgInfo RGBRatioMetric::AnalysisShadow(const cvi* imgShdw, const cvi* imgNonShdw, const cvi* imgShdwMask, 
	OUT vector<double>& gbVs, int histN, int patchRadius)
{
	RGBRatioConverter* converter = new RGBRatioConverter();
	vector<bool> useGain(3, false), useBias(3, false);
	useGain[0] = true; useBias[1] = true; useBias[2] = true;
	return AnalysisShadowBase(imgShdw, imgNonShdw, imgShdwMask, histN, patchRadius, useGain, useBias, -1, -1, gbVs, converter);
}

//---------------------- ShdwAnlysisProc ---------------------------------------

ShdwAnlysisProc::ShdwAnlysisProc(void):m_fileDirs(0)
{
}

ShdwAnlysisProc::~ShdwAnlysisProc(void)
{
}

void ShdwAnlysisProc::AddFileDir(string fileDir)
{
	m_fileDirs.push_back(fileDir);
}

string parGTDir = "GTParam1//";
string shaStaDir = "shadow_analysis//";
int histN = 512;
float resizeRatio = 1.0f/3.0f;
int patchRadius = 7;
float minMaxRatio = 0.49f;

void getParam()
{
	ifstream fin("gtAna.cfg");
	fin>>parGTDir>>shaStaDir>>histN;
	int temp;
	fin>>temp;
	resizeRatio = 1.0f / _f temp;
	fin>>patchRadius;
	fin>>minMaxRatio;
	fin.close();
}

void ShdwAnlysisProc::BunchShdwAnalysis()
{
	RGBGainMetric metric1;
	HLSGainMetric metric2;
	RGBBiasGainMetric metric3;
	LABGainMetric metric4;
	RGBRatioMetric metric5;
	Metric* metric = &metric4;

	getParam();

	ShdwImgInfo sumInfo(histN);
	doFv(k, m_fileDirs)
	{
		cout<<"----------------Handle "<<m_fileDirs[k]<<".---------------\n";
		ShdwImgInfo info = ShdwAnlysis(m_fileDirs[k], histN, resizeRatio, patchRadius, minMaxRatio, metric);
		sumInfo += info;
	}
	cout<<"----------------Overall---------------\n";
// 	sumInfo.computeMinMax(minMaxRatio);
// 	sumInfo.outputMinMax();
// 	cvi* res = VislzDstrbt(sumInfo.biasDstrbt, sumInfo.nPixels);
// 	cvsi("sum_bias.png", res);
// 	cvri(res);
// 	cvi* res2 = VislzDstrbt(sumInfo.gainDstrbt, sumInfo.nPixels);
// 	cvsi("sum_gain.png", res2);
// 	cvri(res2);
// 	cout<<"Image number: "<<sumInfo.nImages<<endl;
// 	cout<<"Avg patch error per image: "<<sumInfo.asmptDist/sumInfo.nImages<<endl;
// 	cout<<"\rDone.";
}

ShdwImgInfo ShdwAnlysisProc::ShdwAnlysis(string fileDir, IN int histN, IN float resizeRatio, IN int patchRadius,
	IN float minMaxRatio, IN Metric* metric)
{
	wGetDirFiles(fileDir + "*_n.png", m_imgNames);
	wMkDir(fileDir + parGTDir);
	//wMkDir(fileDir + shaStaDir);

	ShdwImgInfo sumInfo(histN);
	vector<ShdwImgInfo> info(m_imgNames.size(), histN);

	omp_set_num_threads(1);
#pragma omp parallel for
	doF(i, _i m_imgNames.size())
	{
		cout<<"\rProcessing "<<m_imgNames[i]<<"..";
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		string imgShdwDir = fileDir + imgPrefix + ".png";
		string imgUnShdwDir = fileDir + imgPrefix + "_n.png";
		string imgShdwMaskDir = fileDir + "mask//" + imgPrefix + "s.png";

		cvi* imgShdw = cvlic(imgShdwDir);
		cvi* imgNonShdw = cvlic(imgUnShdwDir);
		cvi* imgShdwMask = cvlig(imgShdwMaskDir);
		cvi* i1 = cvci83(_i (imgShdw->width * resizeRatio), _i (imgShdw->height * resizeRatio));
		cvi* i2 = cvci83(i1), *i3 = cvci81(i1);
		cvResize(imgShdw, i1); cvResize(imgNonShdw, i2); cvResize(imgShdwMask, i3);

		vector<double> gbVs(i1->width * i1->height * 6);

		info[i] = metric->AnalysisShadow(i1, i2, i3, gbVs, histN, patchRadius);
		//info[i] = metric3.AnalysisShadow(imgShdw, imgNonShdw, imgShdwMask, histN);
		info[i].nImages = 1;

		doF(k, 3)
		{
			cvi* v1 = VislzVector(i1->width, i1->height, gbVs, 2*k, 0, 1);
			cvsi(fileDir + parGTDir + imgPrefix + "_" + toStr(k) + "_gain.png", v1);
			cvi* v2 = VislzVector(i1->width, i1->height, gbVs, 2*k+1, -100, 100);
			cvsi(fileDir + parGTDir + imgPrefix + "_" + toStr(k) + "_bias.png", v2);
			cvri(v1); cvri(v2);
			SaveVector(i1->width, i1->height, gbVs, 2*k, fileDir + parGTDir + imgPrefix + "_" + toStr(k) + "_gain.txt");
			SaveVector(i1->width, i1->height, gbVs, 2*k+1, fileDir + parGTDir + imgPrefix + "_" + toStr(k) + "_bais.txt");
		}

// 		cvi* res = VislzDstrbt(info[i].biasDstrbt, info[i].nPixels);
// 		cvsi(fileDir + shaStaDir + imgPrefix + "_bias.png", res);
// 		cvri(res);
// 		cvi* res2 = VislzDstrbt(info[i].gainDstrbt, info[i].nPixels);
// 		cvsi(fileDir + shaStaDir + imgPrefix + "_gain.png", res2);
// 		cvri(res2);

		cvri(imgShdw); cvri(imgNonShdw); cvri(imgShdwMask); cvri(i1); cvri(i2); cvri(i3);
		cout<<m_imgNames[i]<<": "<<info[i].nPixels<<" pixels, avg Patch error = "<<info[i].asmptDist<<"."<<endl;

		info[i].computeMinMax(minMaxRatio);
		info[i].outputMinMax();
	}
	doFv(k, info) sumInfo += info[k];
	//sumInfo.computeMinMax(minMaxRatio);

	cvi* res = VislzDstrbt(sumInfo.biasDstrbt, sumInfo.nPixels);
	cvsi(fileDir + shaStaDir + "sum_bias.png", res);
	cvri(res);
	cvi* res2 = VislzDstrbt(sumInfo.gainDstrbt, sumInfo.nPixels);
	cvsi(fileDir + shaStaDir + "sum_gain.png", res2);
	cvri(res2);
	cout<<"Image number: "<<sumInfo.nImages<<endl;
	cout<<"Avg patch error per image: "<<sumInfo.asmptDist/sumInfo.nImages<<endl;

	//output distribution statistics
	string rangeDir = fileDir + shaStaDir + "range.txt";
	ofstream fout(rangeDir.c_str());
	doF(k, 51)
	{
		sumInfo.computeMinMax(k*0.01f);
		fout<<k*0.01<<":\n";
		doF(k, 3)
		{
			if(sumInfo.useBias[k]) fout<<"Channel "<<k<<" bias: ("<<sumInfo.biasMin[k]<<", "<<sumInfo.biasMax[k]<<").\n";
			if(sumInfo.useGain[k]) fout<<"Channel "<<k<<" gain: ("<<sumInfo.gainMin[k]<<", "<<sumInfo.gainMax[k]<<").\n";
		}
	}
	fout.close();
	pause;

	return sumInfo;
}

void ShdwAnlysisProc::GetStatistics()
{
	string m_fileDir = m_fileDirs[0];
	wGetDirFiles(m_fileDir + "*_n.jpg", m_imgNames);

	float resizeRatio = 0.1f;
	ofstream fout("shadow_analysis//hls.txt");

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

		cvCvtColor(i1, i1, CV_BGR2HLS_FULL);
		cvCvtColor(i2, i2, CV_BGR2HLS_FULL);

		int nPixels = 0;
		//fout<<123<<endl;
		doFcvi(i3, i, j)
		{
			if(cvg20(i3, i, j) < 100) continue;
			cvS v1 = cvg2(i1, i, j), v2 = cvg2(i2, i, j);

			fout<<_i v1.val[0]<<" "<<_i v1.val[1]<<" "<<_i v1.val[2]<<" "
				<<_i v2.val[0]<<" "<<_i v2.val[1]<<" "<<_i v2.val[2]<<endl;
			
			nPixels++;
		}
		cout<<nPixels<<endl;

		cvri(imgShdw); cvri(imgNonShdw); cvri(imgShdwMask); cvri(i1); cvri(i2); cvri(i3);
	}
	fout.close();
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

cvi* ShdwAnlysisProc::VislzFCvi(int width, int height, cvi* src, int ch, double minV, double maxV)
{
	cvi* res = cvci81(width, height);
	doFcvi(res, i, j)
	{
		double v = cvg2(src, i, j).val[ch];
		v = (v - minV) / (maxV - minV) * 255;
		cvs20(res, i, j, v);
	}
	return res;
}

void ShdwAnlysisProc::SaveFCvi(int width, int height, cvi* src, int ch, string fileDir)
{
	ofstream fout(fileDir.c_str());
	fout<<width<<" "<<height<<endl;
	doF(i, height)
	{
		//cout<<"\r"<<i<<" ";
		doF(j, width)
		{ 
			double v = cvg2(src, i, j).val[ch];
			if(fequal(_f v, 0.f)) v = 0;
			fout<<v<<" ";
		}
	}
	fout.close();
}

cvi* ShdwAnlysisProc::VislzVector(int width, int height, vector<double>& vs, int offset, double minV, double maxV)
{
	cvi* res = cvci81(width, height);
	doFcvi(res, i, j)
	{
		double v = vs[(i*width+j)*6+offset];
		v = (v - minV) / (maxV - minV) * 255;
		cvs20(res, i, j, v);
	}
	return res;
}

void ShdwAnlysisProc::SaveVector(int width, int height, vector<double>& vs, int offset, string fileDir)
{
	//cout<<fileDir<<endl;
	ofstream fout(fileDir.c_str());
	fout<<width<<" "<<height<<endl;
	doF(i, height)
	{
		//cout<<"\r"<<i<<" ";
		doF(j, width)
		{ 
			double v = vs[(i*width+j)*6+offset];
			if(fequal(_f v, 0.f)) v = 0;
			fout<<v<<" ";
		}
	}
	fout.close();
}

//---------------------------------------------------------------------------------

void GTCfg::Init()
{
	stepEnabled = true;

	//directory
	srcDir = ".//";
	maskSubDir = "00_Mask//";
	resSubDir = "00_GT//";

	//param
	resizeRatio = 0.5f;
	patchRadius = 7;
}

void GTCfg::InitFromXML(string cfgFile)
{
	tixml::XMLDoc doc;
	if(!doc.Load(cfgFile.c_str())){
		cout<<cfgFile<<" is missing. exit!\n";
		return;
	}

	stepEnabled = (doc.get<int>("gtStep.enabled.@val", 0) == 1);

	//directory
	srcDir = doc.get<string>("srcDir.@val", 0);
	maskSubDir = doc.get<string>("gtStep.maskSubDir.@val", 0);
	resSubDir = doc.get<string>("gtStep.resSubDir.@val", 0);

	//param
	resizeRatio = doc.get<float>("gtStep.resizeRatio.@val", 0);
	patchRadius = doc.get<int>("gtStep.patchRadius.@val", 0);
}

GTCfg ShdwAnlysisProc::cfg;

void ShdwAnlysisProc::LoadCfg(string cfgFile)
{
	cfg.InitFromXML(cfgFile);
}

void ShdwAnlysisProc::GenerateGT()
{
	COutYel("-- GT --\n");
	wMkDir(cfg.srcDir + cfg.resSubDir);
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

	LABGainMetric metric4;
	Metric* metric = &metric4;

	doFv(i, m_imgNames)
	{
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		COutTeal("Processing " + imgPrefix + ".png..\n");

		string imgShdwDir = cfg.srcDir + imgPrefix + ".png";
		string imgUnShdwDir = cfg.srcDir + imgPrefix + "_n.png";
		string imgShdwMaskDir = cfg.srcDir + cfg.maskSubDir + imgPrefix + "s.png";

		cvi* imgShdw = cvlic(imgShdwDir);
		cvi* imgNonShdw = cvlic(imgUnShdwDir);
		cvi* imgShdwMask = cvlig(imgShdwMaskDir);
		if(!imgShdw || !imgNonShdw)
		{
			COutRed("GT image or GT mask is missing!\n");
			cvri(imgShdw); cvri(imgNonShdw); cvri(imgShdwMask); 
			continue;
		}
		if(!imgShdwMask)
		{
			imgShdwMask = cvci81(imgShdw);
			doFcvi(imgShdwMask, i, j) cvs20(imgShdwMask, i, j, 255);
		}

		cvi* i1 = cvci83(_i (imgShdw->width * cfg.resizeRatio), _i (imgShdw->height * cfg.resizeRatio));
		cvi* i2 = cvci83(i1), *i3 = cvci81(i1);
		cvResize(imgShdw, i1); cvResize(imgNonShdw, i2); cvResize(imgShdwMask, i3);

		vector<double> gbVs(i1->width * i1->height * 6);

		metric->AnalysisShadow(i1, i2, i3, gbVs, histN, cfg.patchRadius);

		cvi* param = cvci323(i1);
		doFcvi(param, i, j)
		{
			int idx = 6*(i * param->width + j);
			cvs2(param, i, j, cvs(gbVs[idx], gbVs[idx+3], gbVs[idx+5]));
		}

		ParamLoader::SaveParamToDir(cfg.srcDir + cfg.resSubDir, imgPrefix, param);
		ParamLoader::ShowParamInDir(cfg.srcDir + cfg.resSubDir, imgPrefix, param);

		cvri(imgShdw); cvri(imgNonShdw); cvri(imgShdwMask); cvri(i1); cvri(i2); cvri(i3); cvri(param);
	}
}