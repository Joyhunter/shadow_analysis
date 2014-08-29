#include "StdAfx.h"
#include "recovery_proc.h"
#include "gpm_analysis.h"
#include "shadow_analysis.h"
#include "matrix/matrix.h"
#pragma comment(lib, "sparsemat.lib")

#include "../shadow_removal/mrf_proc.h"
#include "../shadow_removal/decmps_proc.h"

RecoverProc::RecoverProc(void)
{
}

RecoverProc::~RecoverProc(void)
{
}

void RecoverProc::SetFileDir(string fileDir)
{
	m_fileDir = fileDir;
	m_ParamDir = "PredictParam//";
	m_ResultDir = "RecoverResult//";
	m_patchRadius = 3;
}

void RecoverProc::Recover()
{
	wGetDirFiles(m_fileDir + "*_n.png", m_imgNames);
	wMkDir(m_fileDir + m_ResultDir);

	doFv(i, m_imgNames)
	{
		string fileName = m_imgNames[i];
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 6);
		if(imgPrefix != "201") continue;
		cout<<"Handling "<<imgPrefix<<".\n";

		cvi* _srcImg = cvlic(m_fileDir + imgPrefix + ".png");
		cvi* param = GPMAnalysisProc::LoadGTFile(m_fileDir + m_ParamDir, imgPrefix);

		int resizeRatio = 3;
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
			doF(k, 3)
			{
				cvi* temp = ShdwAnlysisProc::VislzFCvi(srcImg->width, srcImg->height, paramResize, k, 
					((k==0)?0:-100), ((k==0)?1:100));
				cvsi(m_fileDir + m_ResultDir + imgPrefix + "_" + toStr(k) + "_o_" + ".png", temp); cvri(temp);
			}
			cvri(paramResize);

			//mrf smoothing and quantify
			cvi* smoothParam = MRFSmooth(srcImg, param);

			//Decmps
			smoothParam = DecmpsSmooth(srcImg, smoothParam); 

			//save smooth param
			doF(k, 3)
			{
				cvi* temp = ShdwAnlysisProc::VislzFCvi(srcImg->width, srcImg->height, smoothParam, k, 
					((k==0)?0:-100), ((k==0)?1:100));
				cvsi(m_fileDir + m_ResultDir + imgPrefix + "_" + toStr(k) + "_s_" + ".png", temp); cvri(temp);
			}

			cvi* recoverRes = ImgRecoverNaive(srcImg, smoothParam);
			//cvi* recoverRes = ImgRecoverPoisson(srcImg, param);
			cvsi(m_fileDir + m_ResultDir + imgPrefix + "__r.png", recoverRes);

			cvri(param); cvri(srcImg); cvri(recoverRes);
			param = smoothParam;
			//resizeRatio /= 2;
			//pause;
		//}
		cvri(_srcImg); cvri(param);
	}
}

cvi* RecoverProc::MRFSmooth(cvi* srcImg, cvi* param)
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
	proc.SolveWithInitial(srcImg, NULL, mask, NULL, 32, resMask);

	doFcvi(mask, i, j)
	{
		auto v = cvg2(paramResize, i, j);
		v.val[0] = cvg20(resMask, i, j) / 255;
		cvs2(paramResize, i, j, v);
	}
	cvri(mask); cvri(resMask);

	return paramResize;
}

cvi* RecoverProc::DecmpsSmooth(cvi* srcImg, cvi* param)
{
	DecmpsProc proc;

	cvi* mask = cvci81(param);
	doFcvi(mask, i, j)
	{
		cvs20(mask, i, j, cvg20(param, i, j)*255);
	}

	cvi* resMask;
	proc.Analysis(srcImg, mask, resMask, 1.0f, 32);

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


cvi* RecoverProc::ImgRecoverNaive(cvi* srcImg, cvi* param)
{
	cvi* res = cvci(srcImg);
	cvZero(res);

	cvi* paramResize = cvci323(srcImg);
	cvResize(param, paramResize);

	cvCvtColor(srcImg, srcImg, CV_BGR2Lab);

	//cvSmooth(paramResize, paramResize, 1, 2*m_patchRadius+1, 2*m_patchRadius+1);
	m_patchRadius = 0;

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
			v.val[1] += v2.val[1];
			v.val[2] += v2.val[2];
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
	cvri(paramResize);

	return res;
}

cvi* RecoverProc::ImgRecoverPoisson(cvi* srcImg, cvi* param)
{
	cvi* res = cvci(srcImg);
	cvZero(res);

	cvi* paramResize = cvci323(srcImg);
	cvResize(param, paramResize);
	cvSmooth(paramResize, paramResize, 1, 2*m_patchRadius+1, 2*m_patchRadius+1);
	m_patchRadius = 0;

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

	//cvCvtColor(res, res, CV_Lab2BGR);
	cvri(paramResize);

	return res;
}