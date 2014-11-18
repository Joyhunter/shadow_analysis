#include "StdAfx.h"
#include "synthesis_proc.h"
#include "matrix/matrix.h"
#include "recovery_proc.h"
#pragma comment(lib, "sparsemat.lib")

SynthesisCfg SynthesisProc::cfg;

SynthesisProc::SynthesisProc(void)
{
	cfg.Init();
}

SynthesisProc::~SynthesisProc(void)
{
}

void SynthesisProc::Test()
{
	string fileDir = "..//dataset//gpmAnalysis//3//synthesis//";
	string imgPrefix = "008";

	cvi* srcImg = cvlic(fileDir + imgPrefix + ".png");
	cvi* holeMask = cvlig(fileDir + imgPrefix + "_hole.png"); cvDilate(holeMask, holeMask, 0, 3);
	cvi* legalMask = cvlig(fileDir + imgPrefix + "_legal.png");
	cvi* guideImg = cvlic(fileDir + imgPrefix + "_guide.png");
	cvi* resImg = cvci(srcImg);

	Synthesis(srcImg, holeMask, resImg, legalMask, guideImg);

	cvsi(fileDir + imgPrefix + "_res.png", resImg);
	cvri(srcImg); cvri(holeMask); cvri(legalMask); cvri(guideImg); cvri(resImg);
}

void SynthesisCfg::Init()
{
	patchSize = 7;

	pymResizeRatio = 0.5f; // 0.5
	pymLevels = 6; // 5 3
	cpltItrN = 20; // 2 30
	cpltItrDecreaseN = 3; // 7 11
	poissonAlpha = 1.0f; // 0.1

	gpmItrN = 6; // 2
	hFlipEnabled = true; // true
	vFlipEnabled = true; // true
	scaleRotateEnabled = false; // false
	gainBiasEnabled = true; // false
	randomSearchEnabled = true;
	randomSearchDistPunish = 1.0f;
	randomSearchScaleFactor = 0.7f;
	randomSearchRadius = 0.2f;
	randomSearchMinRadius = 0.01f;

	guideImgWeight = 100.0f; // 1.0
}

void SynthesisProc::Synthesis(cvi* _srcImg, cvi* _holeMask, cvi* _resImg, cvi* _legalMask, cvi* _guideImg)
{

	float pymResizeRatio = cfg.pymResizeRatio;
	int pymLevels = cfg.pymLevels;
	int cpltItrN = cfg.cpltItrN;
	int gpmItrN = cfg.gpmItrN;
	int patchSize = cfg.patchSize;
	float poissonAlpha = cfg.poissonAlpha;

	int patchOffset = (patchSize - 1) / 2;
	DenseCorrSyn* lastLevelRes = NULL;
	for(int k = 0; k < pymLevels; k++)
	{
		
		if(k == 3)  // 3 4
		{
			//pause;
			cfg.guideImgWeight = 0.1f; //0.1
			cfg.randomSearchDistPunish = 0.3f; //0.3

			cvi* temp = cvlic("syn//1.png"); cvsi("syn//1_prev.png", temp); cvri(temp);
			temp = cvlic("syn//img.png"); cvsi("syn//img_prev.png", temp); cvri(temp);
		}

		//cout<<"Level "<<k<<"...\n";

		//compute new size after resizing
		int resizeRate = _i (1.0f / pymResizeRatio);
		int newW = _srcImg->width, newH = _srcImg->height;
		for(int i = 0; i < pymLevels - 1 - k; i++)
		{
			newW /= resizeRate; newH /= resizeRate;
		}
		int oldW = newW * _i pow(_f resizeRate, pymLevels - 1 - k);
		int oldH = newH * _i pow(_f resizeRate, pymLevels - 1 - k);
		//cout<<newW<<" "<<newH<<endl;

		//doResize and Initialize
		CvRect roi = cvRect(_srcImg->width - oldW, _srcImg->height - oldH, oldW, oldH);
		cvi* img = cvci83(newW, newH);
		cvi* holeMask = cvci81(img); cvSetImageROI(_holeMask, roi); cvResize(_holeMask, holeMask); cvResetImageROI(_holeMask);
		cvi* legalMask = cvci81(img); cvSetImageROI(_legalMask, roi); cvResize(_legalMask, legalMask); cvResetImageROI(_legalMask);
		cvi* guideImg = cvci(img); cvSetImageROI(_guideImg, roi); cvResize(_guideImg, guideImg); cvResetImageROI(_guideImg);
		cvi* resImg = cvci(img);
		InputImageData imageData(img, holeMask, legalMask, guideImg);
		doFcvi(holeMask, i, j)
		{
			if(cvg20(holeMask, i, j) > 0) cvs20(holeMask, i, j, 255);
			if(cvg20(legalMask, i, j) < 255) cvs20(legalMask, i, j, 0);
		}
		//cvErode(legalMask, legalMask, 0, patchOffset + 1);
		if(lastLevelRes == NULL)
		{
			cvSetImageROI(_guideImg, roi);
			cvResize(_guideImg, img);
			cvResetImageROI(_guideImg);
		}
		else
		{
			cvSetImageROI(_srcImg, roi); cvResize(_srcImg, img); cvResetImageROI(_srcImg);
			int boundW = img->width - lastLevelRes->m_width * resizeRate;
			int boundH = img->height - lastLevelRes->m_height * resizeRate;
			lastLevelRes->LevelUp(resizeRate);
			lastLevelRes->AddBound(boundW, 0, boundH, 0);
			lastLevelRes->HandleHoleBoundary(imageData, GPMSynProc::GetRange());
			//lastLevelRes->Identity();
			VoteToCompletion(lastLevelRes, imageData, resImg, patchSize, poissonAlpha);
			//VoteViaLocalCorrection(lastLevelRes, imageData, resImg, patchSize);
			cvCopy(resImg, img);
			lastLevelRes->ShowCorr("syn//1.png");
		}
		cvsi("syn//img.png", img);
		//pause;

		//cvCvtColor(guideImg, guideImg, CV_BGR2Lab);

		//completion iteration
		DenseCorrSyn* dsCor = lastLevelRes;
		for(int l = 0; l < cpltItrN; l++)
		{
			cout<<"\rSyn: Level = "<<k<<", Completion Iteration = "<<l<<"...";
			//if(dsCor != NULL) {delete lastLevelRes; lastLevelRes = NULL;}

			//run gpm
			GPMSynProc proc(1, patchSize, gpmItrN, NULL);
			if(dsCor != NULL) 
			{
				dsCor->UpdatePatchDistance(imageData);
				proc.RunGPMWithInitial(imageData, dsCor);
			}
			else
				dsCor = proc.RunGPM(imageData); //dsCor->ShowCorr("1.png");

			//ofstream fout("1.txt"); dsCor->Save(fout); fout.close();
			//fstream fin("1.txt"); dsCor->Load(fin); fin.close();

			//completion
			VoteToCompletion(dsCor, imageData, resImg, patchSize, poissonAlpha);
			//VoteViaLocalCorrection(dsCor, imageData, resImg, patchSize);
			cvCopy(resImg, img);

			cvsi("syn//img.png", img); //pause;
		}
		//pause;
		//if(lastLevelRes != NULL) {delete lastLevelRes; lastLevelRes = NULL;}
		lastLevelRes = dsCor;

		if(img->width == _resImg->width) cvCopy(img, _resImg);
		cvri(img); cvri(holeMask); cvri(legalMask); cvri(guideImg); cvri(resImg);

		cpltItrN -= cfg.cpltItrDecreaseN;
	}
	cout<<"\rSynthesis complete...\n";
	
	if(lastLevelRes != NULL)
	{
		//cvCopy(lastLevelRes, _resImg);
		delete lastLevelRes; lastLevelRes = NULL;
	}

}

float getUpperIdx(float i, int resizeRatio){return (i+0.5f)*resizeRatio-0.5f;}
float getLowerIdx(float i, int resizeRatio){return (i+0.5f)/resizeRatio-0.5f;}
float getPixelW(int ii, int jj, int i, int j, int patchOffset)
{
	if(patchOffset == 0) return 1;
	return _f gaussian(distEulerL2(_f ii, _f jj, _f i, _f j), 0, _f patchOffset / 2);
}
void getTextureFeature(cvi* src, int i, int j, int patchOffset, cvS& mean, cvS& variance)
{
	mean = cvs(0, 0, 0); variance = cvs(0, 0, 0);
	vector<cvS> patch1; vector<float> weights; float cnt = 0;
	doFs(ii, i-patchOffset, i+patchOffset+1) doFs(jj, j-patchOffset, j+patchOffset+1)
	{
		if(!cvIn(ii, jj, src)) {continue;}
		patch1.push_back(cvg2(src, ii, jj));
		float w = getPixelW(ii, jj, i, j, patchOffset);
		weights.push_back(w); cnt += w;
	}
	doFv(k, patch1){mean += patch1[k] * weights[k];} if(cnt != 0) mean /= cnt;
	doFv(k, patch1){variance += sqr(patch1[k] - mean) * weights[k];}
	if(cnt != 0)doF(k, 3){variance.val[k] = sqrt(variance.val[k]/cnt);}
}
void SynthesisProc::VoteViaLocalCorrection(DenseCorrSyn* dsCor, InputImageData& imageData, cvi* resImg, int patchSize)
{
	int patchOffset = 10;//(patchSize-1) / 2;

	//build pyramid
	int nLevels = 2; float rsRatio = 0.5f; int invRatio = _i(1.f / rsRatio);
	ImgPrmd resPrmd(imageData.guideImg, nLevels, rsRatio), srcPrmd(imageData.src, nLevels, rsRatio);
	//cout<<resImg->width<<endl;

	int resizeRatio = 1;
	doF(k, nLevels+1)
	{
		cvi* levelImg = resPrmd.imgs[k]; //continue;

// 		if(k != nLevels){
// 			resizeRatio *= invRatio;
// 			patchOffset /= invRatio;
// 			continue;
// 		}

		//finer levels
		//1. Compute per-pixel texture feature: mean and sigma
		cvi* meanImg = cvci323(levelImg), *varianceImg = cvci323(levelImg); cvZero(meanImg); cvZero(varianceImg);
		doFcvi(levelImg, i, j)
		{
			//find correspondence position
			float ii = getUpperIdx(_f i, resizeRatio), jj = getUpperIdx(_f j, resizeRatio);
			CorrSyn& corr = dsCor->Get(dsCor->GetCorrIdx(_i ii, _i jj));
			if(cvg20(imageData.holeMask, ii, jj) == 0){corr.x = ii; corr.y = jj;}
			int corI = _i getLowerIdx(corr.x, resizeRatio), corJ = _i getLowerIdx(corr.y, resizeRatio);
			
			//compute mean and variance of non-shadow patch
			cvS avg1 = cvs(0, 0, 0), sigma1 = cvs(0, 0, 0);
			getTextureFeature(levelImg, corI, corJ, patchOffset, avg1, sigma1);
			//avg1 = avg1 * 10 + 128;
			//sigma1 = sigma1 * 10;
			cvs2(meanImg, i, j, avg1); cvs2(varianceImg, i, j, sigma1);
		}
		cvsi("temp3.png", meanImg);
		cvsi("temp4.png", varianceImg);// pause;
		cvsi("temp2.png", levelImg);

		//2. Smooth


		//3. Correct Color, vote back
		cvi* voteSum = cvci323(levelImg), *weight = cvci321(levelImg);
		cvZero(voteSum); cvZero(weight);
		doFcvi(levelImg, i, j)
		{
			//if(i%5!= 0 || j%5!=0) continue;
			float ii = getUpperIdx(_f i, resizeRatio), jj = getUpperIdx(_f j, resizeRatio);
			//if(cvg20(imageData.holeMask, ii, jj) == 0) continue;

			//get texture feature
			cvS mean, variance;
			getTextureFeature(levelImg, i, j, patchOffset, mean, variance);
			cvS mean2 = cvg2(meanImg, i, j), variance2 = cvg2(varianceImg, i, j);

			//vote
			doFs(ii, i-patchOffset, i+patchOffset+1) doFs(jj, j-patchOffset, j+patchOffset+1)
			{
				if(!cvIn(ii, jj, levelImg)) {continue;}

				cvS pixelV = cvg2(levelImg, ii, jj);
				doF(ch, 3)
				{
					if(variance.val[ch] == 0){pixelV.val[ch] = (pixelV.val[ch] - mean.val[ch]) + mean2.val[ch];}
					else
						pixelV.val[ch] = (pixelV.val[ch] - mean.val[ch]) / variance.val[ch] * variance2.val[ch] + mean2.val[ch];
				}
				float pixelW = getPixelW(ii, jj, i, j, patchOffset);

				cvS v = cvg2(voteSum, ii, jj), v2 = cvg2(weight, ii, jj);
				v += pixelV * pixelW; v2 += pixelW;
				cvs2(voteSum, ii, jj, v); cvs2(weight, ii, jj, v2);
			}
		}
		doFcvi(levelImg, i, j)
		{
			cvS v = cvg2(voteSum, i, j), v2 = cvg2(weight, i, j);
			if(v2.val[0] > 1e-3)
			{
				cvs2(levelImg, i, j, v / v2.val[0]); //cout<<v2.val[0]<<" ";
			}
		}
		cvsi("temp1.png", levelImg); pause;

		cvri(meanImg); cvri(varianceImg); cvri(voteSum); cvri(weight);

		resizeRatio *= invRatio;
		patchOffset /= invRatio;
	}

	cvi* corRes = resPrmd.Flatten();
	cvCopy(corRes, resImg);
	resPrmd.Release(); srcPrmd.Release(); cvri(corRes);
}

void SynthesisProc::VoteToCompletion(DenseCorrSyn* dsCor, InputImageData& imageData, cvi* resImg, int patchSize, float gradientAlpha)
{
	int patchOffset = (patchSize-1) / 2;

	cvi* res = cvci323(imageData.src); cvZero(res);
	cvi* gradient1 = cvci323(imageData.src); cvZero(gradient1);
	cvi* gradient2 = cvci323(imageData.src); cvZero(gradient2);
	cvi* cnt = cvci323(imageData.src); cvZero(cnt);
	doFcvi(res, i, j)
	{
		if(cvg20(imageData.holeMask, i, j) == 255) continue; 
		cvs2(res, i, j, cvg2(imageData.guideImg, i, j));
		cvs2(cnt, i, j, cvs(1, 0, 0));
	}
	doFcvi(res, i, j)
	{
		if(cvg20(imageData.holeMask, i, j) == 0) continue;
		//if(!cvIn(i, j, patchOffset, res->height - patchOffset, patchOffset, res->width - patchOffset)) continue;
		CorrSyn& corr = dsCor->Get(dsCor->GetCorrIdx(i, j));
		//if(corr.dist > 100) continue;
		PatchSyn patch;
		PatchDistMetricSyn::GetPatch(imageData, corr.x, corr.y, corr.s, corr.r, corr.hr, corr.vr, (patchSize-1)/2, patch);

		//cvs2(res, i, j, cvg2(imageData.src, corr.x, corr.y));
		//cvs2(cnt, i, j, cvs(1, 1, 1));
		//continue;

		doF(ii, patchSize) doF(jj, patchSize)
		{
			int idx = ii * patchSize + jj;
			int iDel = i + ii - patchOffset, jDel = j + jj - patchOffset;
			if(!cvIn(iDel, jDel, res)) continue;
			if(!cvIn(iDel, jDel, imageData.holeMask)) continue; // || cvg20(imageData.holeMask, iDel, jDel) == 0
			
			cvS v = cvg2(res, iDel, jDel); v += patch.pixels[idx];// * corr.gain + corr.bias;
			cvs2(res, iDel, jDel, v);
			cvS vi = cvg2(cnt, iDel, jDel); vi.val[0]++; cvs2(cnt, iDel, jDel, vi);

			if(ii < patchSize - 1)
			{
				int idx2 = (ii + 1) * patchSize + jj;
				cvS v = cvg2(gradient1, iDel, jDel); v += patch.pixels[idx] - patch.pixels[idx2];
				cvs2(gradient1, iDel, jDel, v);
				cvS vi = cvg2(cnt, iDel, jDel); vi.val[1]++; cvs2(cnt, iDel, jDel, vi);
			}

			if(jj < patchSize - 1)
			{
				int idx2 = ii * patchSize + (jj + 1);
				cvS v = cvg2(gradient2, iDel, jDel); v += patch.pixels[idx] - patch.pixels[idx2];
				cvs2(gradient2, iDel, jDel, v);
				cvS vi = cvg2(cnt, iDel, jDel); vi.val[2]++; cvs2(cnt, iDel, jDel, vi);
			}

		}
	}
	doFcvi(res, i, j)
	{
		//if(cvg20(imageData.holeMask, i, j) == 0) continue;
		cvs2(res, i, j, cvg2(res, i, j) / cvg2(cnt, i, j).val[0]);
		cvs2(gradient1, i, j, cvg2(gradient1, i, j) / cvg2(cnt, i, j).val[1]);
		cvs2(gradient2, i, j, cvg2(gradient2, i, j) / cvg2(cnt, i, j).val[2]);
	}

	cvCopy(imageData.src, resImg);
	//PossionSolve(imageData.src, imageData.holeMask, resImg, res, gradient1, gradient2, gradientAlpha);

	doFcvi(resImg, i, j)
	{
		//if(cvg20(imageData.holeMask, i, j) == 0) continue;
		if(cvg2(cnt, i, j).val[0] == 0) continue;
		cvs2(resImg, i, j, cvg2(res, i, j));
	}

	cvri(res); cvri(cnt); cvri(gradient1); cvri(gradient2);
}

void SynthesisProc::PossionSolve(cvi* srcImg, cvi* holeImg, cvi* resImg, cvi* resConstrain, 
	cvi* gradientConstrainV, cvi* gradientConstrainH, float gradientAlpha)
{
	cvi* holeDilate = cvci(holeImg);
	cvDilate(holeImg, holeDilate);
	cvi* varIdx = cvci321(holeImg); cvZero(varIdx);

	int valN = 0, equN = 0;
	doFcvi(srcImg, i, j)
	{
		if(cvg20(holeDilate, i, j) > 0) 
		{
			valN++; equN++;
			cvs20(varIdx, i, j, valN);
		}
		if(cvg20(holeImg, i, j) > 0)
		{
			if(cvIn(i+1, j, holeImg)) equN++;
			if(cvIn(i, j+1, holeImg)) equN++;
		}
	}
	//cout<<valN<<" "<<equN<<endl;

	sparse::matrix Lmat;
	Lmat.create(equN, valN);
	double* x = new double[valN];
	double* b = new double[equN];
	doF(k, equN) b[k] = 0;

	doF(k, 3)
	{
		int equIdx = 0;
		doFcvi(srcImg, i, j)
		{
			int pIdx = _i cvg20(varIdx, i, j);
			if(pIdx > 0)
			{
				Lmat.add(equIdx, pIdx-1, 1);
				if(cvg20(holeImg, i, j) == 0) b[equIdx++] = cvg2(srcImg, i, j).val[k];
				else b[equIdx++] = cvg2(resConstrain, i, j).val[k];

				if(cvg20(holeImg, i, j) > 0)
				{
					if(cvIn(i+1, j, holeImg))
					{
						int pIdx2 = _i cvg20(varIdx, i+1, j);
						Lmat.add(equIdx, pIdx-1, gradientAlpha); Lmat.add(equIdx, pIdx2-1, -gradientAlpha);
						b[equIdx++] = cvg2(gradientConstrainV, i, j).val[k] * gradientAlpha;
					}
					if(cvIn(i, j+1, holeImg))
					{
						int pIdx2 = _i cvg20(varIdx, i, j+1);
						Lmat.add(equIdx, pIdx-1, gradientAlpha); Lmat.add(equIdx, pIdx2-1, -gradientAlpha);
						b[equIdx++] = cvg2(gradientConstrainH, i, j).val[k] * gradientAlpha;
					}
				}
			}
		}

		Lmat.solve(x, b);

		doFcvi(resImg, i, j)
		{
			int pIdx = _i cvg20(varIdx, i, j);
			if(pIdx > 0)
			{
				double L = x[pIdx-1];
				cvS v = cvg2(resImg, i, j);
				v.val[k] = L;
				cvs2(resImg, i, j, v);
			}
		}
	}

	delete [] x;
	delete [] b;

	cvri(varIdx);
	cvri(holeDilate);
}

//--------------------- Interval & GPMRange ------------------------------------

IntervalSyn::IntervalSyn(float min, float max):min(min), max(max)
{
}

float IntervalSyn::RandValue()
{
	return rand1() * (max - min) + min;
}

//--------------------- Patch & PatchDistMetric ------------------------------------

void PatchDistMetricSyn::GetPatch(InputImageData& src, float x, float y, float s, float r, bool hr, bool vr, 
	int patchOffset, PatchSyn& patch)
{
	int nPixels = sqr(2*patchOffset + 1);

	patch.pixels.resize(nPixels);
	patch.guidePixels.resize(nPixels);

	int idx = 0;
	doFs(pr, -patchOffset, patchOffset+1) doFs(pc, -patchOffset, patchOffset+1)
	{
		int pr2 = hr?(-pr):pr, pc2 = vr?(-pc):pc;
		float preAngle = atan2(_f(pr2), _f(pc2));
		float length = sqrt(sqr(_f(pr2)) + sqr(_f(pc2)));
		float spRowCur = x + length * sin(preAngle - r) * s;
		float spColCur = y + length * cos(preAngle - r) * s;

		if(spRowCur < 0) spRowCur = -spRowCur;
		if(spRowCur > src.src->height - 1) spRowCur = 2 * src.src->height - 2 - spRowCur;
		if(spColCur < 0) spColCur = -spColCur;
		if(spColCur > src.src->width - 1) spColCur = 2 * src.src->width - 2 - spColCur;

		patch.pixels[idx] = cvg2(src.src, spRowCur, spColCur);
		patch.guidePixels[idx] = cvg2(src.guideImg, spRowCur, spColCur);

		idx++;
	}
}

float PatchDistMetricSyn::ComputePatchDist(InputImageData& src, float x, float y, CorrSyn& corr, int patchOffset)
{
	PatchSyn vd, vs;
	GetPatch(src, x, y, 1.f, 0.f, false, false, patchOffset, vd);
	GetPatch(src, corr.x, corr.y, corr.s, corr.r, corr.hr, corr.vr, patchOffset, vs);
	return ComputePatchDistHelper(vd, vs, corr.bias, corr.gain);
}

float PatchDistMetricSyn::ComputePatchDistHelper(PatchSyn& vDst, PatchSyn& vSrc, cvS bias, cvS gain)
{
	float d1 = CptDistDirectWithBiasAndGain(vDst.pixels, vSrc.pixels, cvs(0, 0, 0), cvs(1, 1, 1)); //return d1;
	float d2 = CptDistDirectWithBiasAndGain(vDst.guidePixels, vSrc.guidePixels, bias, gain); //return d2;
	//cvS ratio;
	//float d2 = CptDistAlphaWeight(vDst.guidePixels, vSrc.guidePixels, cvs(-1, -1, -1), cvs(1, 1, 1), false, 3.0f, ratio); //return d2;
	return sqrt(sqr(d1) + sqr(d2)*SynthesisProc::cfg.guideImgWeight);
}

float PatchDistMetricSyn::CptDistDirectWithBiasAndGain(vector<cvS>& vDst, vector<cvS>& vSrc, IN cvS bias, IN cvS gain)
{
	float sum = 0;
	doFv(i, vDst)
	{
		sum += _f cvSDSqr(vDst[i], vSrc[i]*gain+bias);
	}
	return sqrt(sum / vDst.size());
}

float PatchDistMetricSyn::CptDistAlphaWeight(vector<cvS>& vDst, vector<cvS>& vSrc, 
	IN cvS useAlpha, IN cvS weights, IN bool firstIsHue, IN float maxRatio, OUT cvS& alpha)
{
	cvS sumS = cvs(0, 0, 0), sumD = cvs(0, 0, 0);
	doFv(i, vSrc)
	{
		sumS += vSrc[i];
		sumD += vDst[i];
	}
	alpha = cvs(1, 1, 1);
	doF(k, 3)
	{
		float sumS2 = 0, sumD2 = 0;
		doF(i, 3)
		{
			if(useAlpha.val[i] != k) continue;
			sumS2 +=  _f sumS.val[i];
			sumD2 += _f sumD.val[i];
		}
		if(sumD2 > 0) alpha.val[k] = clamp(sumS2 / sumD2, 1.f / maxRatio, maxRatio);
	}

	float sum = 0;
	doFv(i, vDst)
	{
		cvS dstV = vDst[i], srcV = vSrc[i];
		doF(k, 3)
		{
			if(useAlpha.val[k] > 0) dstV.val[k] *= alpha.val[round(useAlpha.val[k])];
		}

		if(firstIsHue)
		{
			float hueChange = fabs(_f (dstV.val[0] - srcV.val[0]));
			if(hueChange > 128) hueChange = 255 - hueChange;
			sum += sqr(hueChange) * _f weights.val[0];
		}
		else
		{
			sum += _f sqr(dstV.val[0] - srcV.val[0]) * _f weights.val[0];
		}

		sum += _f sqr(dstV.val[1] - srcV.val[1]) * _f weights.val[1];
		sum += _f sqr(dstV.val[2] - srcV.val[2]) * _f weights.val[2];
	}

	//doF(i, 3) alphaPenlity *= 2.f - _f gaussian(fabs(_f log(alpha.val[i])), 0.0f, 0.01f);

	return sqrt(sum / vDst.size());
}


//--------------------- Corr & DenseCorr ---------------------------------------

void CorrSyn::Save(ostream& fout)
{
	oswrite(fout, x); oswrite(fout, y);
	oswrite(fout, s); oswrite(fout, r);
	doF(k, 3) oswrite(fout, _f bias.val[k]);
	doF(k, 3) oswrite(fout, _f gain.val[k]);
	oswrite(fout, dist);
}

void CorrSyn::Load(istream& fin)
{
	osread(fin, x); osread(fin, y);
	osread(fin, s); osread(fin, r);
	float temp;
	doF(k, 3)
	{
		osread(fin, temp);
		bias.val[k] = temp;
	}
	doF(k, 3)
	{
		osread(fin, temp);
		gain.val[k] = temp;
	}
	osread(fin, dist);
}

DenseCorrSyn::DenseCorrSyn(int w, int h, int knn, int patchOffset):m_width(w), m_height(h), 
	m_knn(knn), m_patchOffset(patchOffset)
{
	m_values.resize(m_width * m_height * m_knn);
}

CorrSyn DenseCorrSyn::GetRandom(IntervalSyn& hItvl, IntervalSyn& wItvl, GPMSynRange& range, cvi* legalMask)
{
	CorrSyn m_value;
	while(1)
	{
		if(SynthesisProc::cfg.scaleRotateEnabled)
		{
			m_value.x = hItvl.RandValue();
			m_value.y = wItvl.RandValue();
		}
		else
		{
			m_value.x = _f round(hItvl.RandValue());
			m_value.y = _f round(wItvl.RandValue());
		}
		if(!legalMask || cvg20(legalMask, m_value.x, m_value.y) == 255) break;
	}
	m_value.s = range.m_scaleItrl.RandValue();
	m_value.r = range.m_rotateItrl.RandValue();
	doF(p,3) m_value.bias.val[p] = range.m_biasItrl[p].RandValue();
	doF(p,3) m_value.gain.val[p] = range.m_gainItrl[p].RandValue();
	m_value.hr = false; m_value.vr = false;
	if(range.hrEnable) m_value.hr = (rand1() > 0.5f);
	if(range.vrEnable) m_value.vr = (rand1() > 0.5f);
	return m_value;
}

void DenseCorrSyn::RandomInitialize(InputImageData& imgData, GPMSynRange& range)
{
	int idx = -1;

	cvi* legalMask = imgData.legalMask, *holeMask = imgData.holeMask;
	IntervalSyn hItvl(_f m_patchOffset, _f imgData.src->height- 1 - m_patchOffset);
	IntervalSyn wItvl(_f m_patchOffset, _f imgData.src->width - 1 - m_patchOffset);
	doF(i, m_height) doF(j, m_width) doF(k, m_knn)
	{ 
		idx++;
		if(holeMask && cvg20(holeMask, i, j) != 255)
		{
			m_values[idx].dist = -1;
			continue;
		}
		m_values[idx] = GetRandom(hItvl, wItvl, range, legalMask);
	}
	// 	idx = GetCorrIdx(3, 3);
	// 	m_values[idx].x = 3; m_values[idx].y = 3; 
	// 	m_values[idx].s = 1; m_values[idx].r = 0; 
}

void DenseCorrSyn::ShowCorr(string imgStr)
{
	cvi* result = cvci(cvSize(m_width, m_height), 8, 3);
	doFcvi(result, i, j)
	{
		CorrSyn& v = m_values[(i*m_width + j) * m_knn];
		cvs2(result, i, j, cvs(0, v.y / m_width * 255.0, v.x / m_height * 255.0));
	}
	cvsi(imgStr, result);
	cvri(result);
}

void DenseCorrSyn::ShowReflect(string imgStr)
{
	cvi* result = cvci(cvSize(m_width, m_height), 8, 3);
	doFcvi(result, i, j)
	{
		CorrSyn& v = m_values[(i*m_width + j) * m_knn];
		cvs2(result, i, j, cvs(0, v.hr?255:0, v.vr?255:0));
	}
	cvsi(imgStr, result);
	cvri(result);
}

void DenseCorrSyn::ShowCorrDist(string imgStr)
{
	cvi* result = cvci(cvSize(m_width, m_height), 8, 1);
	doFcvi(result, i, j)
	{
		cvs2(result, i, j, cvs(m_values[(i*m_width + j) * m_knn].dist));
	}
	cvsi(imgStr, result);
	cvri(result);
}

void DenseCorrSyn::UpdatePatchDistance(InputImageData& imgData)
{
	int idx = 0;

	doF(i, m_height) doF(j, m_width)
	{
		if(imgData.holeMask && cvg20(imgData.holeMask, i, j) != 255) continue;
		doF(k, m_knn)
		{
			CorrSyn& v = m_values[(i*m_width + j) * m_knn + k];
			v.dist = PatchDistMetricSyn::ComputePatchDist(imgData, _f(i), _f(j), v, m_patchOffset);
			idx++;
		}
		sort(m_values.begin() + (i*m_width + j) * m_knn, 
			m_values.begin() + (i*m_width + j + 1) * m_knn);
	}
}

void DenseCorrSyn::AddCoor(int r, int c, float cx, float cy, float cs, float cr, cvS cBias, cvS cGain, bool hr, bool vr, 
	float cdist)
{
	int sIdx = GetCorrIdx(r, c);
	int eIdx = sIdx + m_knn;
	while(m_values[sIdx].dist < cdist && sIdx < eIdx) sIdx++;
	if(sIdx >= eIdx) return;
	eIdx--;
	while(eIdx > sIdx)
	{
		m_values[eIdx] = m_values[eIdx - 1];
		eIdx--;
	}
	CorrSyn& v = m_values[eIdx];
	v.x = cx; v.y = cy; v.s = cs; v.r = cr; v.hr = hr; v.vr = vr;
	v.bias = cBias; v.gain = cGain;
	v.dist = cdist;
}

void DenseCorrSyn::Save(ostream& fout)
{
	oswrite(fout, m_width); oswrite(fout, m_height);
	oswrite(fout, m_knn); oswrite(fout, m_patchOffset);
	oswrite(fout, _i m_values.size());
	doFv(i, m_values) m_values[i].Save(fout);
}

void DenseCorrSyn::Load(istream& fin)
{
	osread(fin, m_width); osread(fin, m_height);
	osread(fin, m_knn); osread(fin, m_patchOffset);
	int size; osread(fin, size); m_values.resize(size);
	doFv(i, m_values) m_values[i].Load(fin);
}

void DenseCorrSyn::Identity()
{
	doF(i, m_height) doF(j, m_width) doF(k, m_knn)
	{
		int oldIdx = GetCorrIdx(i, j) + k;
		CorrSyn& v = m_values[oldIdx];
		v.x = _f i; v.y = _f j; v.s = 0.5f; v.r = 0.f; v.bias = cvs(0, 0, 0);

		//v.x = (v.x * 0.5f); if(v.x > m_height - 1) v.x -= m_height - 1;
		//v.y = (v.y * 0.5f); if(v.y > m_width - 1) v.y -= m_width - 1;

		v.gain = cvs(1, 1, 1); v.hr = false; v.vr = false;
	}
}

void DenseCorrSyn::LevelUp(int ratio)
{//cout<<m_width<<" "<<m_height<<" "<<ratio<<endl;
	int newWidth = m_width * ratio, newHeight = m_height * ratio;

	vector<CorrSyn> newValues(newWidth * newHeight * m_knn);

	doF(i, m_height) doF(j, m_width) doF(k, m_knn)
	{
		int oldIdx = GetCorrIdx(i, j) + k;
		CorrSyn& oldV = m_values[oldIdx];

		doF(oi, ratio) doF(oj, ratio)
		{
			int idx = ((i * ratio + oi) * newWidth + (j * ratio + oj))*m_knn + k;
			CorrSyn& newV = newValues[idx];

			newV.s = oldV.s;
			newV.r = oldV.r;
			newV.bias = oldV.bias;
			newV.gain = oldV.gain;
			newV.dist = oldV.dist;
			newV.hr = oldV.hr;
			newV.vr = oldV.vr;

			float dx = oi - (_f ratio / 2 - 0.5f);
			float dy = oj - (_f ratio / 2 - 0.5f);
			if(newV.hr) dx = -dx;
			if(newV.vr) dy = -dy;
			float angle = atan2((float)dx, (float)dy);
			float len = pow(dx*dx+dy*dy, 0.5f);
			newV.x = ratio * oldV.x + (_f ratio / 2 - 0.5f) + sin(angle - oldV.r) * oldV.s * len;
			newV.y = ratio * oldV.y + (_f ratio / 2 - 0.5f) + cos(angle - oldV.r) * oldV.s * len;

			//if(oldV.x != 0 && oldV.y != 0)
			//{cout<<oldV.x<<" "<<oldV.y<<" "<<newV.x<<" "<<newV.y<<"   "; pause;}

		}

	}

	m_values = newValues;
	m_width = newWidth;
	m_height = newHeight;
}

void DenseCorrSyn::AddBound(int w1, int w2, int h1, int h2)
{//cout<<m_width<<" "<<m_height<<" "<<w1<<" "<<w2<<" "<<h1<<" "<<h2<<endl;
	int newWidth = m_width + w1 + w2, newHeight = m_height + h1 + h2;
	vector<CorrSyn> newValues(newWidth * newHeight * m_knn);
	doF(i, newHeight) doF(j, newWidth) doF(k, m_knn)
	{
		int oldW = clamp(j - w1, 0, m_width - 1), oldH = clamp(i - h1, 0, m_height - 1);

		int oldIdx = GetCorrIdx(oldH, oldW) + k;
		int idx = (i * newWidth + j) * m_knn + k;
		newValues[idx] = m_values[oldIdx];

		newValues[idx].x += h1;
		newValues[idx].y += w1;
	}
	m_values = newValues;
	m_width = newWidth;
	m_height = newHeight;
	//cout<<m_width<<" "<<m_height<<" "<<w1<<" "<<w2<<" "<<h1<<" "<<h2<<endl;
}

void DenseCorrSyn::HandleHoleBoundary(InputImageData& imgData, GPMSynRange& range)
{
	int idx = -1;
	cvi* legalMask = imgData.legalMask, *holeMask = imgData.holeMask;
	IntervalSyn hItvl(_f m_patchOffset, _f imgData.src->height- 1 - m_patchOffset);
	IntervalSyn wItvl(_f m_patchOffset, _f imgData.src->width - 1 - m_patchOffset);

	doF(i, m_height) doF(j, m_width) doF(k, m_knn)
	{ 
		idx++;
		if(holeMask && cvg20(holeMask, i, j) != 255)
		{
			continue;
		}
		else if(holeMask && cvg20(holeMask, i, j) == 255 && m_values[idx].dist < 0)
		{
			m_values[idx] = GetRandom(hItvl, wItvl, range, legalMask);
			m_values[idx].dist = PatchDistMetricSyn::ComputePatchDist(imgData, _f(i), _f(j), m_values[idx], m_patchOffset);
		}
	}
}

//--------------------- GPMProc ------------------------------------------------

GPMSynRange GPMSynProc::GetRange()
{
	GPMSynRange m_range;

	m_range.hrEnable = SynthesisProc::cfg.hFlipEnabled;//true;
	m_range.vrEnable = SynthesisProc::cfg.vFlipEnabled;//true;

	if(SynthesisProc::cfg.scaleRotateEnabled)
	{
		m_range.setScale(0.67f, 1.5f);
		m_range.setRotate(-1.05f * (float)CV_PI, 1.05f * (float)CV_PI);
	}
	else
	{
		m_range.setScale(1.0f, 1.0f);
		m_range.setRotate(0, 0);
	}

	if(SynthesisProc::cfg.gainBiasEnabled)
	{
		m_range.setGain(cvs(0.9, 0.9, 0.9), cvs(1.1, 1.1, 1.1));
		m_range.setBias(cvs(-5, -5, -5), cvs(5, 5, 5));

// 		m_range.setGain(cvs(0.8, 0.8, 0.8), cvs(1.2, 1.2, 1.2));
// 		m_range.setBias(cvs(-13, -13, -13), cvs(13, 13, 13));
// 
 		m_range.setGain(cvs(0.7, 0.7, 0.7), cvs(1.4, 1.4, 1.4));
 		m_range.setBias(cvs(-20, -20, -20), cvs(20, 20, 20));
	}
	else
	{
		m_range.setGain(cvs(1, 1, 1), cvs(1, 1, 1));
		m_range.setBias(cvs(0, 0, 0), cvs(0, 0, 0));
	}
	return m_range;
}

GPMSynProc::GPMSynProc(int knn, int patchSize, int nItr, GPMSynRange* range): 
m_knn(knn), m_patchSize(patchSize), m_nItr(nItr)
{
	m_patchOffset = (m_patchSize - 1) / 2;

	if(range == NULL)
	{
		m_range = GetRange();
	}
	else
	{
		m_range = *range;
	}

	randInit();
}

GPMSynProc::~GPMSynProc(void)
{
}

DenseCorrSyn* GPMSynProc::RunGPM(InputImageData& src)
{
	int w = src.src->width, h = src.src->height;
	IntervalSyn hItvl(_f m_patchOffset, _f h - 1 - m_patchOffset);
	IntervalSyn wItvl(_f m_patchOffset, _f w - 1 - m_patchOffset);

	int sw = src.src->width, sh = src.src->height;
	IntervalSyn shItvl(_f m_patchOffset, _f sh - 1 - m_patchOffset);
	IntervalSyn swItvl(_f m_patchOffset, _f sw - 1 - m_patchOffset);

	//random initialize
	DenseCorrSyn* dsCor = new DenseCorrSyn(w, h, m_knn, m_patchOffset);
	dsCor->RandomInitialize(src, m_range);
	dsCor->UpdatePatchDistance(src);
	//dsCor->ShowCorr("1.png");
	//dsCor->ShowReflect("2.png");
	//dsCor->ShowCorrDist("initDist.png");

	RunGPMWithInitial(src, dsCor);

	return dsCor;
}

void GPMSynProc::RunGPMWithInitial(InputImageData& src, DenseCorrSyn* dsCor)
{
// 	int w = src.src->width, h = src.src->height;
// 	IntervalSyn hItvl(_f m_patchOffset, _f h - 1 - m_patchOffset);
// 	IntervalSyn wItvl(_f m_patchOffset, _f w - 1 - m_patchOffset);
	int sw = src.src->width, sh = src.src->height;
	IntervalSyn shItvl(_f m_patchOffset, _f sh - 1 - m_patchOffset);
	IntervalSyn swItvl(_f m_patchOffset, _f sw - 1 - m_patchOffset);

	int w = src.src->width, h = src.src->height;
	IntervalSyn hItvl(0, _f h - 1);
	IntervalSyn wItvl(0, _f w - 1);

	//iteration
	doF(kItr, m_nItr)
	{
		int nCores = 1;
		doF(k, nCores)
		{
			float hStart = hItvl.min, hEnd = hItvl.max;
			float wStart = wItvl.min, wEnd = wItvl.max;
			if(m_nItr % 4 == 0 || m_nItr % 4 == 1)
			{
				float hD = (hEnd - hStart) / nCores;
				hStart = hD * k + hStart;
				hEnd = hStart + hD;
			}
			else
			{
				float wD = (wEnd - wStart) / nCores;
				wStart = wD * k + wStart;
				wEnd = wStart + wD;
			}

			int propSum = 0, ranSum = 0;
			doFs(_x, _i hStart, _i hEnd) doFs(_y, _i wStart, _i wEnd)
			{
				//if(_y == wStart) cout<<"\rInteration "<<kItr<<": "<<(_x-hStart)*100/(hEnd-hStart)<<"% done.";
				int x = _x, y = _y;
				if(kItr % 2 == 1){
					x = _i hItvl.min + _i hItvl.max - x;
					y = _i wItvl.min + _i wItvl.max - y;
				}
				//if((kItr % 4 == 1) || (kItr % 4 == 2))
				if(src.holeMask && cvg20(src.holeMask, x, y) != 255) continue;
				PatchSyn dstPatch;
				PatchDistMetricSyn::GetPatch(src, _f x, _f y, 1.f, 0.f, false, false, m_patchOffset, dstPatch);
				propSum += Propagate(src, x, y, (kItr % 2 == 1), *dsCor, wItvl, hItvl, 
					swItvl, shItvl, dstPatch);
				if(SynthesisProc::cfg.randomSearchEnabled)
					ranSum += RandomSearch(src, x, y, *dsCor, swItvl, shItvl, dstPatch);
			}
			//cout<<"\nPropSum = "<<propSum<<", RanSum = "<<ranSum<<".\n";
		}
		dsCor->ShowCorr("syn//1.png");
		dsCor->ShowReflect("syn//2.png");
	}
	//cout<<"\rGPM complete.\n";
}

int GPMSynProc::Propagate(InputImageData& src, int x, int y, bool direction, 
	DenseCorrSyn& dsCor, IntervalSyn wItvl, IntervalSyn hItvl, 
	IntervalSyn swItvl, IntervalSyn shItvl, PatchSyn& dstPatch)
{
	int offset[4][2] = {{1, 0}, {0, 1}, {0, -1}, {-1, 0}};
	int offsetIdxStart = 0;
	int offsetIdxEnd = 2;
	if(!direction){
		offsetIdxStart = 2;
		offsetIdxEnd = 4;
	}

	int propSum = 0;

	float distThres = dsCor.GetDistThres(x, y);
	doFs(i, offsetIdxStart, offsetIdxEnd)
	{
		int dx = offset[i][0], dy = offset[i][1];

		if(!cvIn(x + dx, y + dy, _i hItvl.min, _i hItvl.max, 
			_i wItvl.min, _i wItvl.max))
			continue;

		if(src.holeMask && cvg20(src.holeMask, x + dx, y + dy) != 255) continue;

		int corrOffIdx = dsCor.GetCorrIdx(x + dx, y + dy);
		doF(k, m_knn)
		{
			CorrSyn& v = dsCor.Get(corrOffIdx + k); //v.dist *= 2;
			float angle = atan2(-(float)dx, -(float)dy);
			float newX = clamp(v.x + sin(angle - v.r) * v.s, shItvl.min, shItvl.max);
			float newY = clamp(v.y + cos(angle - v.r) * v.s, swItvl.min, swItvl.max);
			
			if(v.hr) newX = clamp(v.x - sin(angle - v.r) * v.s, shItvl.min, shItvl.max);
			if(v.vr) newY = clamp(v.y - cos(angle - v.r) * v.s, swItvl.min, swItvl.max);

			//cout<<v.x<<" "<<v.y<<" "<<newX<<" "<<newY<<"   ";

			if(src.legalMask && cvg20(src.legalMask, newX, newY) != 255) continue;

			PatchSyn srcPatch;
			PatchDistMetricSyn::GetPatch(src, newX, newY, v.s, v.r, v.hr, v.vr, m_patchOffset, srcPatch);
			float patchDist = PatchDistMetricSyn::ComputePatchDistHelper(dstPatch, srcPatch, v.bias, v.gain);
			//patchDist *= SynthesisProc::cfg.randomSearchDistPunish;

			if(patchDist * SynthesisProc::cfg.randomSearchDistPunish  < distThres){
				propSum++;
				dsCor.AddCoor(x, y, newX, newY, v.s, v.r, v.bias, v.gain, v.hr, v.vr, patchDist);
				distThres = dsCor.GetDistThres(x, y);
			}
		}
	}
	return propSum;
}

int GPMSynProc::RandomSearch(InputImageData& src, int x, int y, 
	DenseCorrSyn& dsCor, IntervalSyn wItvl, IntervalSyn hItvl, PatchSyn& dstPatch)
{
	float distThres = dsCor.GetDistThres(x, y);
	int idx = dsCor.GetCorrIdx(x, y);

	int ranSum = 0;

	vector<CorrSyn> corrs(m_knn);
	doF(k, m_knn) corrs[k] = dsCor.Get(idx + k);

	float scaleFactor = SynthesisProc::cfg.randomSearchScaleFactor;
	float searchRadius = SynthesisProc::cfg.randomSearchRadius;
	doF(k, m_knn)
	{
		CorrSyn& v = corrs[k];
		float hSpace = hItvl.max - hItvl.min;
		float wSpace = wItvl.max - wItvl.min;
		float sSpace = m_range.m_scaleItrl.max - m_range.m_scaleItrl.min;
		float rSpace = m_range.m_rotateItrl.max - m_range.m_rotateItrl.min;

		hSpace = hSpace * searchRadius; wSpace = wSpace * searchRadius; sSpace *= searchRadius; rSpace *= searchRadius;

		vector<float> biasSpace(m_range.m_biasItrl.size());
		doFv(k, m_range.m_biasItrl) biasSpace[k] = (m_range.m_biasItrl[k].max - m_range.m_biasItrl[k].min) * searchRadius;
		vector<float> gainSpace(m_range.m_gainItrl.size());
		doFv(k, m_range.m_gainItrl) gainSpace[k] = (m_range.m_gainItrl[k].max - m_range.m_gainItrl[k].min) * searchRadius;

		while(hSpace > (hItvl.max - hItvl.min) * SynthesisProc::cfg.randomSearchMinRadius)
		{
			//cout<<"\r"<<hSpace<<" "<<wSpace;
			float newx = clamp(hSpace * (2 * rand1() - 1) + v.x, hItvl.min, hItvl.max);
			float newy = clamp(wSpace * (2 * rand1() - 1) + v.y, wItvl.min, wItvl.max);
			
			if(!SynthesisProc::cfg.scaleRotateEnabled)
			{
				newx = _f round(newx); newy = _f round(newy);
			}

			if(!src.legalMask || cvg20(src.legalMask, newx, newy) == 255)
			{
				float news = clamp(sSpace * (2 * rand1() - 1) + v.s, m_range.m_scaleItrl.min, m_range.m_scaleItrl.max);
				float newr = clamp(rSpace * (2 * rand1() - 1) + v.r, m_range.m_rotateItrl.min, m_range.m_rotateItrl.max);
				cvS newBias, newGain;
				doF(k, 3)
				{
					newBias.val[k] = clamp(biasSpace[k] * (2 * rand1() - 1) + _f v.bias.val[k], 
						m_range.m_biasItrl[k].min, m_range.m_biasItrl[k].max);
					newGain.val[k] = clamp(gainSpace[k] * (2 * rand1() - 1) + _f v.gain.val[k], 
						m_range.m_gainItrl[k].min, m_range.m_gainItrl[k].max);
				}
				bool newHr = (rand1() > 0.5f);
				bool newVr = (rand1() > 0.5f);

				PatchSyn srcPatch;
				PatchDistMetricSyn::GetPatch(src, newx, newy, news, newr, newHr, newVr, m_patchOffset, srcPatch);
				float patchDist = PatchDistMetricSyn::ComputePatchDistHelper(dstPatch, srcPatch, newBias, newGain);

				if(patchDist < distThres * SynthesisProc::cfg.randomSearchDistPunish){
					ranSum++;
					dsCor.AddCoor(x, y, newx, newy, news, newr, newBias, newGain, newHr, newVr, patchDist);
					distThres = dsCor.GetDistThres(x, y);
				}
			}

			hSpace *= scaleFactor; wSpace *= scaleFactor;
			sSpace *= scaleFactor; rSpace *= scaleFactor;
			doFv(k, m_range.m_biasItrl) biasSpace[k] *= scaleFactor;
			doFv(k, m_range.m_gainItrl) gainSpace[k] *= scaleFactor;
		}
	}
	return ranSum;
}