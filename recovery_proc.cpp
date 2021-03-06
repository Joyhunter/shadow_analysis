#include "StdAfx.h"
#include "recovery_proc.h"
#include "param_loader.h"
#include "shadow_analysis.h"
#include "matrix/matrix.h"
#pragma comment(lib, "sparsemat.lib")

#include "../shadow_removal/mrf_proc.h"
#include "../shadow_removal/decmps_proc.h"

//----------------------------- Recovery Configure ----------------------------

void RecoveryCfg::Init()
{
	stepEnabled = true;

	//directory
	srcDir = ".//";
	paramSubdir = "02_VoteDetection//";
	resSubdir = "03_VoteRecovery//";

	//focused files
	focusedPrefix.clear();
	focusedPrefix.push_back("012");

	//step 1: compute naive res
	patchsize_naiveRecovery = 1;

	//step 2: synthesis
	//mask generation
	shdwDegreeThres_maskGenerate = 0.8f;
	shdwBoundDilateTimes_maskGenerate = 0.01f;
	//synthesis
	SynthesisProc::cfg.Init();
	gtAsGuidance = false;

	//step 3: local color correction
	patchRadius_localColorCorrection = 10;
	correctionStepRatio_localColorCorrection = 0.5f;
	//pyramid
	pyramidExtraLevels_pyramid = 3;
	pyramidResizeRatio_pyramid = 0.5f;

	naiveResImgPostfix = "_naive_res.png";
	holeMaskImgPostfix = "_syn_hole.png";
	legalMaskImgPostfix = "_syn_legal.png";
	synResImgPostfix = "_syn_res.png";
	corResImgPostfix = "_cor_res.png";

	//origial
	patchRadius_old = 3;
}

void RecoveryCfg::InitFromXML(string cfgFile)
{
	tixml::XMLDoc doc;
	if(!doc.Load(cfgFile.c_str())){
		cout<<cfgFile<<" is missing. exit!\n";
		return;
	}

	stepEnabled = (doc.get<int>("recoveryStep.enabled.@val", 0) == 1);

	//directory
	srcDir = doc.get<string>("srcDir.@val", 0);
	paramSubdir = doc.get<string>("detectionStep.subDir.@val", 0); 
	resSubdir = doc.get<string>("recoveryStep.subDir.@val", 0);

	//focused files
	int num = doc.get<int>("focusedImgs.@num", 0);
	focusedPrefix.clear();
	doF(k, num)	focusedPrefix.push_back(doc.get<string>("focusedImgs.focuedPrefix.@val", k));

	//step 1: compute naive res
	patchsize_naiveRecovery = doc.get<int>("recoveryStep.naiveRecovery.patchSize.@val", 0);

	//step 2: synthesis
	//mask generation
	shdwDegreeThres_maskGenerate = doc.get<float>("recoveryStep.maskGenerate.shdwDegreeThres.@val", 0);
	shdwBoundDilateTimes_maskGenerate = doc.get<float>("recoveryStep.maskGenerate.shdwBoundDilateTimes.@val", 0);
	//synthesis
	gtAsGuidance = (doc.get<int>("recoveryStep.synthesis.gtAsGuidance.@val", 0) == 1);

	//step 3: local color correction
	patchRadius_localColorCorrection = doc.get<int>("recoveryStep.localColorCorrection.patchRadius.@val", 0);
	correctionStepRatio_localColorCorrection = doc.get<float>("recoveryStep.localColorCorrection.correctionStepRatio.@val", 0);
	//pyramid
	pyramidExtraLevels_pyramid = doc.get<int>("recoveryStep.localColorCorrection.pyramid.extraLevels.@val", 0);
	pyramidResizeRatio_pyramid = doc.get<float>("recoveryStep.localColorCorrection.pyramid.resizeRatio.@val", 0);

	naiveResImgPostfix = doc.get<string>("recoveryStep.naiveRecovery.resPostfix.@val", 0);
	holeMaskImgPostfix = doc.get<string>("recoveryStep.maskGenerate.holeMaskImgPostfix.@val", 0);
	legalMaskImgPostfix = doc.get<string>("recoveryStep.maskGenerate.legalMaskImgPostfix.@val", 0);
	synResImgPostfix = doc.get<string>("recoveryStep.synthesis.resPostfix.@val", 0);
	corResImgPostfix = doc.get<string>("recoveryStep.localColorCorrection.resPostfix.@val", 0);

	//origial
	patchRadius_old = doc.get<int>("recoveryStep.patchRadius_old.@val", 0);

	SynthesisProc::cfg.InitFromXML(cfgFile);

// 	cout<<patchsize_naiveRecovery<<" "<<shdwDegreeThres_maskGenerate<<" "<<shdwBoundDilateTimes_maskGenerate<<" "<<patchRadius_localColorCorrection
// 		<<" "<<correctionStepRatio_localColorCorrection<<" "<<pyramidExtraLevels_pyramid<<" "<<pyramidResizeRatio_pyramid<<" "
// 		<<patchRadius_old;
}

//----------------------------- Set Configure File -----------------------------

RecoveryCfg RecoverProc::cfg;

RecoverProc::RecoverProc(void)
{
	cfg.Init();
}

RecoverProc::~RecoverProc(void)
{
}

void RecoverProc::LoadCfg(string cfgFile)
{
	cfg.InitFromXML(cfgFile);
}

//----------------------------- Recovery ----------------------------------------

void RecoverProc::Recover3()
{
	COutYel("-- Recovery --\n");
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

		//load source image and param prediction
		cvi* srcImg = cvlic(cfg.srcDir + imgPrefix + ".png");
		cvi* freeImg = cvlic(cfg.srcDir + imgPrefix + "_n.png");
		cvi* param = ParamLoader::LoadParamFromDir(cfg.srcDir + cfg.paramSubdir, imgPrefix);

// 		cvi* param = cvci323(srcImg);
// 		cvResize(_param, param);

		//Step 1: Get naive recovery result
		COutput("Start step 1: get naive result..", CC_DARKGREEN);
		cvi* naiveRes = ImgRecoverNaive(srcImg, param, cfg.patchsize_naiveRecovery);
		if(cfg.naiveResImgPostfix != "") cvsi(cfg.srcDir + cfg.resSubdir + imgPrefix + cfg.naiveResImgPostfix, naiveRes);
		COutGreen("\r                                   \rStep 1: Naive result got.\n");

		//cvi* localRes = cvlic(m_fileDir + m_ResultDir + imgPrefix + "_syn_res.png");
		//Step 2: Get Patch Synthesis result
		COutput("Generating shadow hole and legal mask..", CC_DARKGREEN);
		cvi* holeMask = cvci81(srcImg), *legalMask = cvci81(srcImg);
		GenerateMaskFromParam(param, holeMask, legalMask, cfg.shdwDegreeThres_maskGenerate, cfg.shdwBoundDilateTimes_maskGenerate);
		if(cfg.holeMaskImgPostfix != "") cvsi(cfg.srcDir + cfg.resSubdir + imgPrefix + cfg.holeMaskImgPostfix, holeMask);
		if(cfg.legalMaskImgPostfix != "") cvsi(cfg.srcDir + cfg.resSubdir + imgPrefix + cfg.legalMaskImgPostfix, legalMask);
		COutGreen("\rStart step 2: Hole and legal mask got, start synthesis...\n");
		SynthesisProc sProc;
		cvi* synRes = cvci(srcImg);
		if(cfg.gtAsGuidance && freeImg) sProc.Synthesis(srcImg, holeMask, synRes, legalMask, freeImg);
		else
		{
			if(cfg.gtAsGuidance && !freeImg) COutRed("Free image load failed, use naive res as guidance.\n");
			sProc.Synthesis(srcImg, holeMask, synRes, legalMask, naiveRes);
		}
		if(cfg.synResImgPostfix != "") cvsi(cfg.srcDir + cfg.resSubdir + imgPrefix + cfg.synResImgPostfix, synRes);
		COutGreen("\rStep 2: Synthesis result got.\n");
// 		cvi* synRes = cvlic(m_fileDir + m_ResultDir + imgPrefix + "_syn_res.png");
// 		cvi* naiveRes = cvlic(m_fileDir + m_ResultDir + imgPrefix + "_naive_res.png");
// 		cvi* holeMask = NULL; cvi* param = NULL; cvi* legalMask = NULL;

		//Step 3: Local color correction
		COutput("Start step 3: local color correction..", CC_DARKGREEN);
		cvi* corRes = LocalColorCorrection(naiveRes, synRes, holeMask);
		if(cfg.corResImgPostfix != "") cvsi(cfg.srcDir + cfg.resSubdir + imgPrefix + cfg.corResImgPostfix, corRes);
		COutGreen("\r                                       \rStep 3: Correction result got.\n");

		cvri(srcImg); cvri(param); cvri(naiveRes); cvri(synRes); cvri(freeImg); //cvri(localRes);
		cvri(holeMask); cvri(legalMask); cvri(corRes);
	}
}

cvi* RecoverProc::ImgRecoverNaive(cvi* _srcImg, cvi* param, int patchSize)
{
	cvi* res = cvci(_srcImg);
	cvZero(res);

	cvi* paramResize = cvci323(_srcImg);
	cvResize(param, paramResize);

	cvi* srcImg = cvci(_srcImg);
	cvCvtColor(srcImg, srcImg, CV_BGR2Lab);

	//cvSmooth(paramResize, paramResize, 1, 2*m_patchRadius+1, 2*m_patchRadius+1);
	int m_patchRadius = (patchSize - 1)/2;

	vector<cvS> resVec(srcImg->width * srcImg->height, cvs(0, 0, 0));
	vector<int> counts(srcImg->width * srcImg->height, 0);

	doFcvi(srcImg, i, j)
	{
		cvS v2 = cvg2(paramResize, i, j);

		doFs(ii, -m_patchRadius, m_patchRadius+1) doFs(jj, -m_patchRadius, m_patchRadius+1)
		{
			if(!cvIn(i+ii, j+jj, srcImg)) continue;
			cvS v = cvg2(srcImg, i+ii, j+jj);
			v.val[0] /= v2.val[0];

			//v.val[1] = v.val[0]; v.val[2] = v.val[0];
			v.val[1] += v2.val[1]; v.val[2] += v2.val[2];

			int idx = (i+ii)*srcImg->width + j+jj;
			resVec[idx] += v;
			counts[idx]++;
		}
	}

	doFcvi(srcImg, i, j)
	{
		int idx = (i)*srcImg->width + j;
		cvs2(res, i, j, resVec[idx] / counts[idx]);
	}

	cvCvtColor(res, res, CV_Lab2BGR);
	cvri(paramResize); cvri(srcImg);

	return res;
}

void RecoverProc::GenerateMaskFromParam(IN cvi* param, OUT cvi* holeMask, OUT cvi* legalMask, float shadowThres, float boundW)
{
	cvZero(holeMask); cvZero(legalMask);

	doFcvi(param, i, j)
	{
		if(cvg20(param, i, j) < shadowThres)
			cvs20(holeMask, i, j, 255);
	}

	int w = param->width;
	cvDilate(holeMask, holeMask, 0, _i(_f w*boundW));
	doFcvi(legalMask, i, j)
	{
		cvs20(legalMask, i, j, 255 - cvg20(holeMask, i, j));
	}
}

cvi* RecoverProc::GetCorCpltPixels(cvi* _naiveRes, cvi* _synRes)
{
	cvi* mask = cvci83(_naiveRes); cvZero(mask);

	cvi* naiveRes = cvci(_naiveRes); cvCvtColor(naiveRes, naiveRes, CV_BGR2Lab);
	cvi* synRes = cvci(_synRes); cvCvtColor(synRes, synRes, CV_BGR2Lab);

	int patchOffset = 2;

	doFcvi(mask, i, j)
	{
		cvS avg1 = cvs(0, 0, 0), avg2 = avg1;
		int n = 0;
		doFs(ii, -patchOffset, patchOffset + 1) doFs(jj, -patchOffset, patchOffset + 1)
		{
			if(!cvIn(i+ii, j+jj, mask)) continue;
			n++;
			avg1 += cvg2(naiveRes, i+ii, j+jj); avg2 += cvg2(synRes, i+ii, j+jj);
		}
		avg1 /= n; avg2 /= n;

		//decide
		double lGain = (avg1.val[0] < 20)?1:10; if(avg2.val[0] != 0) lGain = avg1.val[0] / avg2.val[0];
		double aBias = fabs(avg1.val[1] - avg2.val[1]), bBias = fabs(avg1.val[2] - avg2.val[2]);

		cvs2(mask, i, j, cvs(lGain*100, aBias, bBias));
		//if(lGain > 0.8f && lGain < 1.2f && aBias < 3 && bBias < 3) cvs20(mask, i, j, 255);
	}
	cvri(naiveRes); cvri(synRes);
	
	return mask;
}

cvi* RecoverProc::LocalColorCorrection(cvi* _naiveRes, cvi* _synRes, cvi* _holeMask)
{	
	int patchRadius = cfg.patchRadius_localColorCorrection;

	//get correct completion pixels
	cvi* cpltMask = GetCorCpltPixels(_naiveRes, _synRes);
	//cvsi("4.png", cpltMask);

	//build pyramid
	int nLevels = cfg.pyramidExtraLevels_pyramid;
	float rsRatio = cfg.pyramidResizeRatio_pyramid;
	ImgPrmd naivePrmd(_naiveRes, nLevels, rsRatio), synPrmd(_synRes, nLevels, rsRatio);

	doF(k, nLevels+1)
	{
		cvi* naiveRes = naivePrmd.imgs[k], *synRes = synPrmd.imgs[k];
		cvi* holeMask = NULL; if(_holeMask){holeMask = cvci81(naiveRes); cvResize(_holeMask, holeMask);}

		cvi* corRes = LocalColorCorrectionSingleLevel(naiveRes, synRes, holeMask, patchRadius);
		//cvsi("1.png", naiveRes); cvsi("2.png", synRes); cvsi("3.png", corRes); //pause;

		cvCopy(corRes, naiveRes); cvri(corRes);		
		cvri(holeMask);

		patchRadius = max2(_i (_f patchRadius * rsRatio), 1); 
	}

	cvi* corRes = naivePrmd.Flatten();
	naivePrmd.Release(); synPrmd.Release();

	cvri(cpltMask);

	return corRes;
}

cvi* RecoverProc::LocalColorCorrectionSingleLevel(cvi* naiveRes, cvi* synRes, cvi* holeMask, int patchRadius)
{
	cvi* corRes = cvci(naiveRes); //cvZero(corRes);

	int step = max2(_i(patchRadius * cfg.correctionStepRatio_localColorCorrection), 1);

	cvi* voteSum = cvci323(corRes), *weight = cvci321(corRes);
	cvZero(voteSum); cvZero(weight);
	doFcvi(naiveRes, i, j)
	{
		if(holeMask && cvg20(holeMask, i, j) == 0) continue;
		if(i % step != 0 || j % step != 0) continue;

		//get patches
		vector<cvS> patch1, patch2; vector<float> weights; float cnt = 0;
		doFs(ii, i-patchRadius, i+patchRadius+1) doFs(jj, j-patchRadius, j+patchRadius+1)
		{
			if(!cvIn(ii, jj, corRes)) {continue;}
			patch1.push_back(cvg2(naiveRes, ii, jj)); patch2.push_back(cvg2(synRes, ii, jj));
			float w = _f gaussian(distEulerL2(_f ii, _f jj, _f i, _f j), 0, _f patchRadius / 2);
			weights.push_back(w); cnt += w;
		}

		//get avg and sigma
		cvS avg1 = cvs(0, 0, 0), avg2 = cvs(0, 0, 0), sigma1 = cvs(0, 0, 0), sigma2 = cvs(0, 0, 0);
		doFv(k, patch1){ avg1 += patch1[k] * weights[k]; avg2 += patch2[k] * weights[k]; }
		avg1 /= cnt; avg2 /= cnt;
		doFv(k, patch1){ sigma1 += sqr(patch1[k] - avg1) * weights[k]; sigma2 += sqr(patch2[k] - avg2) * weights[k];}
		doF(k, 3){ sigma1.val[k] = sqrt(sigma1.val[k]/cnt); sigma2.val[k] = sqrt(sigma2.val[k]/cnt);}
		//cout<<avg1.val[0]<<" "<<avg2.val[0]<<" "<<sigma1.val[0]<<" "<<sigma2.val[0]<<"  "; pause;
		doFv(k, patch1) doF(ch, 3)
		{
			if(sigma1.val[ch] == 0){patch1[k].val[ch] = (patch1[k].val[ch] - avg1.val[ch]) + avg2.val[ch];}
			else
				patch1[k].val[ch] = (patch1[k].val[ch] - avg1.val[ch]) / sigma1.val[ch] * sigma2.val[ch] + avg2.val[ch];
		}

		//set back
		int idx = 0;
		doFs(ii, i-patchRadius, i+patchRadius+1) doFs(jj, j-patchRadius, j+patchRadius+1)
		{
			if(!cvIn(ii, jj, corRes)) {continue;}
			cvS v = cvg2(voteSum, ii, jj), v2 = cvg2(weight, ii, jj);
			v += patch1[idx] * weights[idx]; v2 += weights[idx];
			//v += patch1[idx]; v2 += 1;
			cvs2(voteSum, ii, jj, v); cvs2(weight, ii, jj, v2);
			idx++;
		}
	}

	//output vote Result
	doFcvi(corRes, i, j)
	{
		cvS v = cvg2(voteSum, i, j), v2 = cvg2(weight, i, j);
		if(v2.val[0] > 1e-3)
		{
			cvs2(corRes, i, j, v / v2.val[0]); //cout<<v2.val[0]<<" ";
		}
	}//cvsi(corRes);

	cvri(voteSum); cvri(weight);
	return corRes;
}


//----------------------------- Image Pyramid ----------------------------------------

ImgPrmd::ImgPrmd(cvi* img, int levels, float ratio)
{
	imgs.clear();
	cvi* src = cvci323(img);
	doFcvi(src, i, j) cvs2(src, i, j, cvg2(img, i, j));

	doF(k, levels)
	{
		cvi* ri = cvci323(_i(src->width*ratio), _i(src->height*ratio));
		cvResize(src, ri);
		cvi* rir = cvci323(src);
		cvResize(ri, rir);

		cvi* res = cvci323(src);
		doFcvi(res, i, j)
			cvs2(res, i, j, cvg2(src, i, j) - cvg2(rir, i, j));

		imgs.push_back(res);
		cvri(src); cvri(rir);
		src = ri;
		//cvsi(res); pause;
	}

	imgs.push_back(src);
}
cvi* ImgPrmd::Flatten()
{
	for (size_t i = imgs.size() - 1; i > 0; i--)
	{
		cvi* low = imgs[i], *high = imgs[i-1];
		cvi* low2 = cvci323(high);
		cvResize(low, low2);
		doFcvi(low2, i, j)
			cvs2(high, i, j, cvg2(high, i, j) + cvg2(low2, i, j));
		cvri(low2);
	}

	cvi* res = cvci83(imgs[0]);
	doFcvi(res, i, j) cvs2(res, i, j, cvg2(imgs[0], i, j));
	return res;
}
void ImgPrmd::Release()
{
	doFv(k, imgs) cvri(imgs[k]);
}

cvi* RecoverProc::GetBestParamLapPymd(IN cvi* src, IN cvi* mask, IN cvi* cpImg, OUT cvS& param)
{
	int nLevels = 3; float rsRatio = 0.5f;
	ImgPrmd srcPrmd(src, nLevels, rsRatio), cpPrmd(cpImg, nLevels, rsRatio);

	cvi* mask32 = cvci321(mask);
	doFcvi(mask, i, j) cvs20(mask32, i, j, 1.0 - cvg20(mask, i, j) / 255);

	float minP = 0.4f, maxP = 1.5f; 
	int sampleN = 50;

	doF(k, nLevels)
	{
		//if(k != nLevels) continue;
		cvi* s = srcPrmd.imgs[k], *c = cpPrmd.imgs[k];
		cvi* m = cvci321(s); cvResize(mask32, m);

		//cvsi(s); pause; cvsi(c); pause; cvsi(m); pause;
		//cvsi("__1.png", s); cvsi("__2.png", c); cvsi("__3.png", m); pause;

		if(k != nLevels)
		{
			float athres = 0.2f;
			doF(kk, 3)
			{
				float avgS = 0, avgC = 0, sigS = 0, sigC = 0, nP = 0;
				doFcvi(s, i, j)
				{
					float w = _f cvg20(m, i, j); if(w < athres) w = 0; else w = 1;
					//if(w < athres) continue;
					avgS += _f cvg2(s, i, j).val[kk] * w;
					avgC += _f cvg2(c, i, j).val[kk] * w; nP+= w;
				}
				avgS /= nP; avgC /= nP;
				doFcvi(s, i, j)
				{
					float w = _f cvg20(m, i, j); if(w < athres) w = 0; else w = 1;
					//if(w < athres) continue;
					sigS += pow(_f cvg2(s, i, j).val[kk] - avgS, 2.0f) * w;
					sigC += pow(_f cvg2(c, i, j).val[kk] - avgC, 2.0f) * w;
				}
				sigS = pow(sigS / nP, 0.5f);
				sigC = pow(sigC / nP, 0.5f); 
				//cout<<nP<<" "<<avgS<<" "<<avgC<<" "<<sigS<<" "<<sigC<<endl;
				doFcvi(s, i, j)
				{
					if(cvg20(m, i, j) < athres) continue;
					auto vs = cvg2(s, i, j);
					double r = (vs.val[kk] - avgS) / sigS;
					vs.val[kk] = (r * sigC) + avgC;
					cvs2(s, i, j, vs);
				}
			}
		}
		else
		{
			doF(kk, 3)
			{
				float minError = 1e10, minRatio = 0;
				doF(si, sampleN)
				{ 
					float ratio = (maxP - minP) / (sampleN - 1) * si + minP;
					float avgError = 0;
					doFcvi(s, i, j)
					{
						float v = _f cvg2(s, i, j).val[kk];
						float w = _f cvg20(m, i, j);
						avgError += (_f cvg2(c, i, j).val[kk] - (v / (1.0f - w*ratio))) * w;
					}
					avgError = fabs(avgError);
					avgError /= s->width * s->height;


					if(avgError < minError)
					{
						minError = avgError; minRatio = ratio;
					}
				}
				doFcvi(s, i, j)
				{
					auto v = cvg2(s, i, j);
					float w = _f cvg20(m, i, j);
					v.val[kk] /= (1.0 - w*minRatio);
					cvs2(s, i, j, v);
				}
				//cout<<minRatio<<" ";
			}
		}
		//cvsi("__1.png", s); cvsi("__2.png", c); cvsi("__3.png", m);
		//pause;


		cvri(m);
		//pause;
	}

	cvi* res = srcPrmd.Flatten();
	cvsi(res);
	
	srcPrmd.Release(); cpPrmd.Release();

	return res;
}

cvi* RecoverProc::GetBestParam(IN cvi* src, IN cvi* mask, IN cvi* cpImg, OUT cvS& param)
{
	cvi* res = cvci(src);

	float minP = 0.2f, maxP = 2.0f; 
	int sampleN = 50;

	//whether a pixel is legal
	cvi* legalMap = cvci81(src); cvZero(legalMap);
	doFcvi(legalMap, i, j)
	{
		cvS v1 = cvg2(src, i, j), v2 = cvg2(cpImg, i, j);
		v1 = v1/(v2+0.01);
		if(v1.val[0] > 0.4 && v1.val[0] < 2.0 && 
			v1.val[1] > 0.4 && v1.val[1] < 2.0 &&
			v1.val[2] > 0.4 && v1.val[2] < 2.0)
			cvs20(legalMap, i, j, 255);
	}
	cvsi(legalMap);

	cvi* t1 = cvci81(src);
	doF(k, 3)
	{
		float minSig = 1e10, minRatio = 0;
		
// 		float Lsrc = 0, Lcp = 0, wSum = 0;
// 		doFcvi(src, i, j)
// 		{
// 			float v = _f cvg2(src, i, j).val[k];
// 			float v2 = _f cvg2(cpImg, i, j).val[k];
// 			float v3 = 1.0 - _f cvg20(mask, i, j) / 255.0f;
// 
// 			Lsrc += v * v3;
// 			Lcp += v2 * v3;
// 			wSum += (1.0 - v3) * v3;
// 		}
// 		float factor = Lsrc / Lcp ;


		doF(si, sampleN)
		{ 
			float ratio = (maxP - minP) / (sampleN - 1) * si + minP;
			doFcvi(src, i, j)
			{
				float v = _f cvg2(src, i, j).val[k];
				float v2 = _f cvg20(mask, i, j) / 255.0f;
				v2 = 1.0f - (1.0f - v2) * ratio;
				cvs20(t1, i, j, v/v2);
			}

			//compute sigma
			float sigma = 0;
// 			float avg = 0;
// 			doFcvi(t1, i, j) avg += _f cvg20(t1, i, j);
// 			avg /= t1->width * t1->height;
// 			doFcvi(t1, i, j) sigma += pow(_f cvg20(t1, i, j) - avg, 2.0f);
			doFcvi(t1, i, j)
			{
				//if(cvg20(legalMap, i, j) == 0) continue;
				float v3 = 1.0f - _f cvg20(mask, i, j) / 255.0f;
				sigma += (_f cvg2(cpImg, i, j).val[k] - _f cvg20(t1, i, j)) * v3;
			}
			sigma = fabs(sigma);
 			//sigma /= t1->width * t1->height;


			if(sigma < minSig)
			{
				minSig = sigma; minRatio = ratio;
			}
		}
		cout<<minRatio<<" "; //pause;

		doFcvi(res, i, j)
		{
			auto v = cvg2(res, i, j);
			float v2 = _f cvg20(mask, i, j) / 255.0f;
			v2 = 1.0f - (1.0f - v2) * minRatio;
			v.val[k] /= v2;
			cvs2(res, i, j, v);
		}
		cvsic(res, k);
		//pause;

	}

	return res;
}

//-----------------------------------------------------------------------------------

//currently not use
cvi* RecoverProc::BoundarySmooth(cvi* srcImg, cvi* smoothParam, cvi* &shadowMask, float Lthres, float r1, float r2)
{
	//get boundary Mask
	cvi* mask = cvci81(smoothParam); cvZero(mask);
	doFcvi(mask, i, j)
	{
		if(cvg20(smoothParam, i, j)*255.f < Lthres) cvs20(mask, i, j, 255); 
	}
	cvi* mask2 = cvci(mask);
	shadowMask = cvci(mask);
	cvDilate(shadowMask, shadowMask, 0, _i r2);

	//get All levels pixel
	int levelsN = _i r1 + _i r2;
	vector<vector<pair<int, int> > > allLevelpixels(levelsN);
	cvi* maskBak = cvci(mask);
	doF(k, _i r1)
	{
		cvErode(mask, mask2, 0, 1);
		doFcvi(mask, i, j)
		{
			if(cvg20(mask, i, j) != cvg20(mask2, i, j)) allLevelpixels[_i r1-1-k].push_back(make_pair(i, j));
		}
		cvCopy(mask2, mask);
	}
	cvi* boundMask = cvci(mask);
	cvCopy(maskBak, mask); cvri(maskBak);
	doF(k, _i r2)
	{
		cvDilate(mask, mask2, 0, 1);
		doFcvi(mask, i, j)
		{
			if(cvg20(mask, i, j) != cvg20(mask2, i, j)) allLevelpixels[k+_i r1].push_back(make_pair(i, j));
		}
		cvCopy(mask2, mask);
	}
	doFcvi(mask2, i, j)
	{
		if(cvg20(mask2, i, j) == 255 && cvg20(boundMask, i, j) == 0) cvs20(boundMask, i, j, 255);
		else cvs20(boundMask, i, j, 0);
	} return boundMask;

	//all level avg Lumi
	cvCvtColor(srcImg, srcImg, CV_BGR2Lab);
	vector<cvS> avgEachLOri(levelsN), avgEachLRev(levelsN);
	doF(k, levelsN)
	{
		cvS avgOri = cvs(0, 0, 0), avgRev = cvs(0, 0, 0);
		for(auto p = allLevelpixels[k].begin(); p != allLevelpixels[k].end(); p++)
		{
			auto v = cvg2(srcImg, p->first, p->second);
			auto v2 = cvg2(smoothParam, p->first, p->second);
			avgOri += v;
			v.val[0] /= v2.val[0]; v.val[1] += v2.val[1]; v.val[2] += v2.val[2];
			avgRev += v;
		}
		avgEachLOri[k] = avgOri / allLevelpixels[k].size();
		avgEachLRev[k] = avgRev / allLevelpixels[k].size();
		//cout<<avgEachLOri[k].val[0]<<" "<<avgEachLRev[k].val[0]<<endl;
	}
	cvS step = (avgEachLRev[levelsN - 1] - avgEachLRev[0]) / (levelsN - 1);
	doF(k, levelsN)
	{
		cvS oldParam = avgEachLRev[k] - avgEachLOri[k];
		oldParam.val[0] = avgEachLOri[k].val[0] / avgEachLRev[k].val[0];
		avgEachLRev[k] = step * k + avgEachLRev[0];
		cvS newParam = avgEachLRev[k] - avgEachLOri[k];
		newParam.val[0] = avgEachLOri[k].val[0] / avgEachLRev[k].val[0];
		//cout<<oldParam.val[0]<<" "<<newParam.val[0]<<endl;

		for(auto p = allLevelpixels[k].begin(); p != allLevelpixels[k].end(); p++)
		{
			cvS v = cvg2(smoothParam, p->first, p->second);
			v.val[0] = v.val[0] / oldParam.val[0] * newParam.val[0];
			v.val[1] = v.val[1] - oldParam.val[1] + newParam.val[1];
			v.val[2] = v.val[2] - oldParam.val[2] + newParam.val[2];
			cvs2(smoothParam, p->first, p->second, v);
		}
	}


	cvCvtColor(srcImg, srcImg, CV_Lab2BGR);

	cvri(mask); cvri(mask2);
	return boundMask;
}

cvi* RecoverProc::VislzBound(cvi* src, cvi* boundMask)
{
	cvi* res = cvci(src);
	doFcvi(res, i, j)
	{
		if(cvg20(boundMask, i, j) == 0)
			cvs2(res, i, j, cvg2(src, i, j) / 5);
	}
	return res;
}

cvi* RecoverProc::ImgRecoverPoisson(cvi* srcImg, cvi* param)
{
	cvi* res = cvci(srcImg);
	cvZero(res);

	cvi* paramResize = cvci323(srcImg);
	cvResize(param, paramResize);
	cvSmooth(paramResize, paramResize, 1, 2*cfg.patchRadius_old+1, 2*cfg.patchRadius_old+1);
	cfg.patchRadius_old = 0;

	cvCvtColor(srcImg, srcImg, CV_BGR2Lab);

	//solve poisson
	int pixelSize = srcImg->width * srcImg->height;
	double al1 = 0;
	sparse::matrix Lmat;	
	Lmat.create(pixelSize*5, pixelSize);
	double* x = new double[pixelSize];
	double* b = new double[5*pixelSize];
	doF(k, 5*pixelSize) b[k] = 0;

	int idx = 0;
	doFcvi(srcImg, i, j)
	{
		cvS v = cvg2(srcImg, i, j), v2 = cvg2(paramResize, i, j);
		if(v2.val[0] > 0.99)
		{
			Lmat.add(5*idx, idx, 0.0001);
			b[5*idx] = v.val[0] * 0.0001;
		}
		else
		{
			Lmat.add(5*idx, idx, al1);
			b[5*idx] = v.val[0] / v2.val[0] * al1;

			int off = 1;
			doFs(ii, -1, 2) doFs(jj, -1, 2)
			{
				if(abs(ii) + abs(jj) != 1) continue;
				if(!cvIn(i+ii, j+jj, srcImg)) continue;
				cvS v3 = cvg2(srcImg, i+ii, j+jj);
				int idx2 = (i+ii)*srcImg->width + j+jj;
				Lmat.add(5*idx+off, idx, 1); Lmat.add(5*idx+off, idx2, -1);
				b[5*idx+off] = (v.val[0] - v3.val[0]) / v2.val[0];
				off++;
			}
		}

		idx++;
	}

	Lmat.solve(x, b);

	doFcvi(res, i, j)
	{
		double L = x[i*res->width + j];
		cvs2(res, i, j, cvs(L, L, L));
	}

	delete [] x;
	delete [] b;

	cvCvtColor(res, res, CV_Lab2BGR);
	cvri(paramResize);

	return res;
}

void RecoverProc::PoissonSmooth(cvi* srcImg, cvi* mask)
{
	cvCvtColor(srcImg, srcImg, CV_BGR2Lab);

	//solve poisson
	int pixelSize = srcImg->width * srcImg->height;
	double al1 = 0;
	sparse::matrix Lmat;	
	Lmat.create(pixelSize*5, pixelSize);
	double* x = new double[pixelSize];
	double* b = new double[5*pixelSize];
	doF(k, 5*pixelSize) b[k] = 0;

	int idx = 0;
	doFcvi(srcImg, i, j)
	{
		if(cvg20(mask, i, j) == 0)
		{
			bool isEdge = false;
			doFs(ii, -1, 2) doFs(jj, -1, 2)
			{
				if(cvIn(i+ii, j+jj, srcImg) && cvg20(mask, i+ii, j+jj) > 100)
					isEdge = true;
			}
			if(isEdge)
			{
				Lmat.add(5*idx, idx, 1);
				b[5*idx] = cvg20(srcImg, i, j);
			}
		}
		else
		{
			int off = 1;
			doFs(ii, -1, 2) doFs(jj, -1, 2)
			{
				if(abs(ii) + abs(jj) != 1) continue;
				if(!cvIn(i+ii, j+jj, srcImg)) continue;
				int idx2 = (i+ii)*srcImg->width + j+jj;
				Lmat.add(5*idx+off, idx, 1); Lmat.add(5*idx+off, idx2, -1);
				off++;
			}
		}
		idx++;
	}

	Lmat.solve(x, b);

	doFcvi(srcImg, i, j)
	{
		if(cvg20(mask, i, j) == 0) continue;
		auto v = cvg2(srcImg, i, j);
		v.val[0] = x[i*srcImg->width + j];
		cvs2(srcImg, i, j, v);
	}

	delete [] x;
	delete [] b;

	cvCvtColor(srcImg, srcImg, CV_Lab2BGR);
}

//abandoned
void RecoverProc::SetFileDir(string fileDir)
{
// 	m_fileDir = fileDir;
// 	m_ParamDir = "02_Detection//"; // VoteParam PredictParam
// 	m_ResultDir = "03_Recovery//";
// 
// 	m_ParamDir = "02_VoteDetection//"; // VoteParam PredictParam
// 	m_ResultDir = "03_VoteRecovery//";
// 
// 	m_patchRadius = 3;
}

void RecoverProc::GetScript()
{
	ofstream fout("test.bat");
	doFv(i, m_imgNames)
	{
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		//fout<<"PMRefImpl.exe rgb8 "<<imgPrefix<<"__.png "<<imgPrefix<<"res.png -hole_mask "<<imgPrefix<<"__m.png"
		//<<" -init_guess_image "<<imgPrefix<<"__r4.png -em_iters_min 50 -em_iters_final 50 -min_coarse_size 75\n";
		fout<<"PMRefImpl.exe rgb8 "<<imgPrefix<<"__r3.png "<<imgPrefix<<"__c3.png -hole_mask ..//"<<imgPrefix<<"__r3.png"
			<<" -init_guess_image "<<imgPrefix<<"__r1.png -em_iters_min 50 -em_iters_final 50 -min_coarse_size 75\n";
	}
	fout.close();
}



/*

void RecoverProc::Recover()
{
	wGetDirFiles(m_fileDir + "*_n.png", m_imgNames);
	wMkDir(m_fileDir + m_ResultDir);

	//GetScript(); return;

	doFv(i, m_imgNames)
	{
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		if(imgPrefix != "004") continue;
		cout<<"Handling "<<imgPrefix<<".\n";

		cvi* _srcImg = cvlic(m_fileDir + imgPrefix + ".png");
		cvi* param = GPMAnalysisProc::LoadGTFile(m_fileDir + m_ParamDir, imgPrefix);

		int resizeRatio = 1;
		//while(1)
		//{
		cvi* srcImg = cvci83(_srcImg->width / resizeRatio, _srcImg->height / resizeRatio);
		cvResize(_srcImg, srcImg);
		//cvi* temp = cvci(srcImg);
		//cvSmooth(temp, srcImg, CV_BILATERAL, 20, 20, 200, 200); cvri(temp);

		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__.png", srcImg);

		//save init param
		cvi* paramResize = cvci323(srcImg);
		cvResize(param, paramResize);
		// 			doF(k, 3)
		// 			{
		// 				cvi* temp = ShdwAnlysisProc::VislzFCvi(srcImg->width, srcImg->height, paramResize, k, 
		// 					((k==0)?0:-100), ((k==0)?1:100));
		// 				cvsi(m_fileDir + m_ResultDir + imgPrefix + "_" + toStr(k) + "_o_" + ".png", temp); cvri(temp);
		// 			}
		cvri(paramResize);

		//mrf smoothing and quantify
		cvi* smoothParam = MRFSmooth(srcImg, param, 16);
		smoothParam = DecmpsSmooth(srcImg, smoothParam, 128);

		//before boundary smooth save
		// 			doF(k, 3)
		// 			{
		// 				cvi* temp = ShdwAnlysisProc::VislzFCvi(srcImg->width, srcImg->height, smoothParam, k, 
		// 					((k==0)?0:-100), ((k==0)?1:100));
		// 				cvsi(m_fileDir + m_ResultDir + imgPrefix + "_" + toStr(k) + "_s1_" + ".png", temp); cvri(temp);
		// 			}
		cvi* recoverRes = ImgRecoverNaive(srcImg, smoothParam);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__r1.png", recoverRes);

		//boundary Smooth
		cvi* shadowBinaryMask;
		float lthres = 185;
		if(imgPrefix == "114") lthres = 210;
		if(imgPrefix == "116") lthres = 230;
		if(imgPrefix == "109") lthres = 220;
		if(imgPrefix == "102") lthres = 235;
		if(imgPrefix == "100") lthres = 210;
		//if(imgPrefix == "201") lthres = 230;
		cvi* boundMask = BoundarySmooth(srcImg, smoothParam, shadowBinaryMask, 
			lthres, srcImg->width * 0.02f, srcImg->width * 0.02f); //185
		cvi* boundShow = VislzBound(recoverRes, boundMask);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__r2.png", boundShow);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__m.png", shadowBinaryMask);
		cvri(boundShow); cvri(shadowBinaryMask);

		//after boundary smooth save smooth param
		doF(k, 3)
		{
			cvi* temp = ShdwAnlysisProc::VislzFCvi(srcImg->width, srcImg->height, smoothParam, k, 
				((k==0)?0:-100), ((k==0)?1:100));
			cvsi(m_fileDir + m_ResultDir + imgPrefix + "_" + toStr(k) + "_s2_" + ".png", temp); cvri(temp);
		}
		cvi* recoverRes2 = ImgRecoverNaive(srcImg, smoothParam); 
		//cvi* recoverRes = ImgRecoverPoisson(srcImg, param);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__r3.png", boundMask);

		//Smooth boundary lightness
		PoissonSmooth(recoverRes2, boundMask);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__r4.png", recoverRes2);

		cvri(param); cvri(srcImg); cvri(recoverRes); cvri(recoverRes2); cvri(boundMask);
		param = smoothParam;
		//resizeRatio /= 2;
		//pause;
		//}
		cvri(_srcImg); cvri(param);
	}
}

void RecoverProc::Recover2()
{
	wGetDirFiles(m_fileDir + "*_n.png", m_imgNames);
	wMkDir(m_fileDir + m_ResultDir);

	//GetScript(); return;

	doFv(i, m_imgNames)
	{
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		if(imgPrefix != "316") continue;
		cout<<"Handling "<<imgPrefix<<".\n";

		cvi* srcImg = cvlic(m_fileDir + imgPrefix + ".png");
		cvi* _maskImg = cvlig(m_fileDir + m_ResultDir + imgPrefix + "_0_s2_.png");
		cvi* maskImg = cvci81(srcImg); cvResize(_maskImg, maskImg); cvri(_maskImg);
		//cvi* maskImg = cvlig(m_fileDir + "GTParam3//" + imgPrefix + "_0_gain.png");
		//cvSmooth(maskImg, maskImg, 2, srcImg->width / 40, srcImg->width / 40);
		//cvi* cpImg = cvlic(m_fileDir + m_ResultDir + "pm_cp_res//NEng_nRNG_Win_x64_Rel_8_" + imgPrefix + "res.png");
		cvi* _cpImg = cvlic(m_fileDir + imgPrefix + "_n.png"); cvi* cpImg = cvci83(srcImg); cvResize(_cpImg, cpImg); cvri(_cpImg);

		//match source with cp
		//laplacian pyramid, only use correct completion result
		cvS gain = cvs(0, 0, 0);
		cvi* recoverImg = GetBestParam(srcImg, maskImg, cpImg, gain);
		cvsi(m_fileDir + m_ResultDir + "results//" + imgPrefix + "__r1.png", recoverImg);
		recoverImg = GetBestParamLapPymd(recoverImg, maskImg, cpImg, gain); //pause;
		cvsi(m_fileDir + m_ResultDir + "results//" + imgPrefix + "__r2.png", recoverImg);

		//get boundary
		float thres = 185;
		float ratio1 = 0.02f, ratio2 = 0.02f;
		if(imgPrefix == "004"){thres = 150; ratio1 = 0.01f; ratio2 = 0.04f;}
		if(imgPrefix == "012"){ratio1 = 0.01f; ratio2 = 0.04f;}
		if(imgPrefix == "015"){ratio1 = 0.03f; ratio2 = 0.03f;}
		// 		if(imgPrefix == "114"){thres = 210;}
		// 		if(imgPrefix == "102"){thres = 242;}
		// 		if(imgPrefix == "116") thres = 230;
		// 		if(imgPrefix == "109") thres = 220;
		// 		if(imgPrefix == "100") thres = 210;
		// 		if(imgPrefix == "307"){thres = 245; ratio1 = 0.05f; ratio2 = 0.05f;}
		cvi* shadowBinaryMask;
		cvi* boundMask = BoundarySmooth(srcImg, maskImg, shadowBinaryMask, 
			thres*255, srcImg->width * ratio1, srcImg->width * ratio2);
		cvi* boundShow = VislzBound(srcImg, boundMask);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__m.png", shadowBinaryMask);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__b_show.png", boundShow);
		cvsi(m_fileDir + m_ResultDir + imgPrefix + "__b.png", boundMask);
		//break;

		//poisson boundary
		//PoissonSmooth(recoverImg, boundMask);
		doFcvi(recoverImg, i, j) if(cvg20(boundMask, i, j) > 100) cvs2(recoverImg, i, j, cvg2(srcImg, i, j));
		cvri(boundShow); cvri(boundMask); cvri(shadowBinaryMask);

		cvsi(m_fileDir + m_ResultDir + "results//" + imgPrefix + "__r3.png", recoverImg);


		cvri(srcImg); cvri(maskImg); cvri(recoverImg);


	}
}

cvi* RecoverProc::MRFSmooth(cvi* srcImg, cvi* param, int nLabels)
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

cvi* RecoverProc::DecmpsSmooth(cvi* srcImg, cvi* param, int nLabels)
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

*/