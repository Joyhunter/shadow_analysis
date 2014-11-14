#include "StdAfx.h"
#include "gpm_analysis.h"
#include "shadow_analysis.h"

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

void GPMAnalysisProc::ReadParam(string cfgFile)
{
	ifstream fin(cfgFile.c_str());
	float v1, v2, v3, v4, v5, v6;
	fin>>v1>>v2;
	param.range.setScale(v1, v2);
	fin>>v1>>v2;
	param.range.setRotate(v1 * (float)CV_PI, v2 * (float)CV_PI);
	fin>>v1>>v2>>v3>>v4>>v5>>v6;
	param.range.setGain(cvs(v1, v2, v3), cvs(v4, v5, v6)); //0.92
	fin>>v1>>v2>>v3>>v4>>v5>>v6;
	param.range.setBias(cvs(v1, v2, v3), cvs(v4, v5, v6)); //21.4
	fin>>param.downRatio;
	fin>>param.multiLevelRatio;
	fin>>param.extraLevelN;
	fin>>param.patchSize;
	fin>>param.distThres;
	fin>>rangeDir;
	fin.close();
}

string GPMAnalysisProc::GetCorrFileDir(string imgName)
{
	return "corr//" + rangeDir + imgName.substr(0, imgName.size()-4)
		+ "_" + toStr(param.downRatio) + "_" + toStr(param.multiLevelRatio) + "_" + toStr(param.extraLevelN)
		+ "_" + toStr(param.patchSize) + ".corr";
}

//----------------------------- interface 1: Run Grid GPM for all images -------------------

void GPMAnalysisProc::RunGPMForAllImages()
{
	cout<<"Reading image names...";
	wGetDirFiles(m_fileDir + "*.png", m_imgNames);

	cout<<"Reading params...";
	ReadParam();

	cout<<"Begin Running gpm for all images..."<<endl;

	wMkDir(m_fileDir + "corr//" + rangeDir);

	doFv(i, m_imgNames)
	{
		if(m_imgNames[i][m_imgNames[i].size()-5] == 'n') continue;
		cout<<"\r  Runing GPM for image "<<m_imgNames[i]<<"...\n";

		ImgContainer img(m_fileDir + m_imgNames[i], param.downRatio, param.colorMode);

		int gridSize = _i(img.src()->width * 0.28f);
		int gridOffset = _i(img.src()->width * 0.25f);

		string coorName = GetCorrFileDir(m_imgNames[i]);
		string corrFileName = m_fileDir + coorName;

		LmnIvrtPatchDistMetric metric2;
		GridGPMProc proc(&metric2, gridSize, gridOffset, 1, param.patchSize, 20, &(param.range));

		proc.RunGridGPMMultiScale(img, corrFileName, param.multiLevelRatio, param.extraLevelN);
		//if(img.srcR() == NULL) img.GenerateResizedImg(1);
		//proc.ShowGPMResUI(img, corrFileName, param.distThres);

		cout<<"\r  GPM of "<<m_imgNames[i]<<" complete. File saved to "<<coorName<<".\n";
	}
}

//----------------------------------------------------------------------------------------------











void GPMAnalysisProc::ShdwAnlysis()
{
	wGetDirFiles(m_fileDir + "*.png", m_imgNames);

	//0.01: (0.29, 0.94). (-5.47, 8.47). (-6.47461, 21.416)
	//0.05: (0.391602, 0.827148). (-3.48633, 6.47461). (-3.48633, 18.4277).
	//0.1:  (0.43457, 0.780273). (-2.49023, 4.48242). (-2.49023, 17.4316).
	//0.15: (0.463867, 0.743164). (-1.49414, 4.48242). (-0.498047, 16.4355).
	//0.2: (0.493164, 0.717773). (-1.49414, 3.48633). (0.498047, 15.4395).
	//0.25: (0.510742, 0.698242). (-0.498047, 3.48633). (4.48242, 14.4434).
	//0.3: (0.522461, 0.682617). (-0.498047, 2.49023). (6.47461, 13.4473).
	//0.4: (0.549805, 0.649414). (0.498047, 2.49023). (7.4707, 11.4551).
	//0.49: (0.594727, 0.608398). (1.49414, 1.49414). (9.46289, 9.46289).

	//range 1(0.01): 0.26 0.87  -6.48 8.47   -5.5 19.5
	//range 2(0.1): 0.43 0.78  -2.49 4.48   -2.49 17.43
	//range 3(0.15): 0.46 0.74  -1.49 4.49   -0.49 16.44
	//range 4(0.2): 0.49 0.72  -1.49 3.49  0.49 15.44  set range wrong..
	//range 5(0.25): 0.51 0.70  -0.49 3.49  4.48 14.44
	//range 6(0.3): 0.52 0.68  -0.49 2.49  6.47 13.45
	//range 7(0.4): 0.54 0.65  0.49 2.49  7.47 11.46
	//range 8(0.49): 0.59 0.61  1.49 1.50  9.46 9.47 
	//range 9(0.5): 0.60 0.60  1.50 1.50  9.46 9.46

	param.range.setScale(0.67f, 1.5f);
	param.range.setRotate(-1.05f * (float)CV_PI, 1.05f * (float)CV_PI);
	param.range.setGain(cvs(0.6, 1, 1), cvs(0.6, 1, 1)); //0.92
	param.range.setBias(cvs(0, 1.5, 9.46), cvs(0, 1.5, 9.46)); //21.4
	param.downRatio = 6; param.multiLevelRatio = 2; param.extraLevelN = 0;
	param.patchSize = 7;
	param.colorMode = CV_BGR2Lab;
	param.distThres = 30;

	rangeDir = "range4//";
	parVoteDir = "VoteParam//";
	parGTDir = "GTParam6//";

	//ReadParam();

// 	if(1)
// 	{
// 		//RunGPMForAllImages(); return;
// 
// 		ofstream fout("sta.txt");
// 		doF(i, 128)
// 		{
// 			cout<<"\r"<<i;
// 			param.distThres = _f 2 * i;
// 			GetStatistics(fout);
// 		}
// 		fout.close();
// 	}

	test(); return;

// 	VoteInitMask(); return;
// 	ofstream fout2("voteError.txt");
// 	doF(i, 40)
// 	{
// 		cout<<"\r"<<i*2;
// 		param.distThres = _f 2 * i;
// 		VoteInitMask(); 
// 		fout2<<GetVoteError()<<endl;
// 	}
// 	fout2.close();
// 	return;
}

void GPMAnalysisProc::test()
{
	ifstream fin("vote.cfg");
	fin>>rangeDir;
	fin>>parVoteDir;
	fin>>parGTDir;
	fin>>param.downRatio;
	fin>>param.multiLevelRatio;
	fin>>param.extraLevelN;
	fin>>param.distThres;
	fin.close();
	wMkDir(m_fileDir + parVoteDir);
	//VoteInitMask(); GetVoteError(); pause; return;
	ofstream fout2("voteError.txt");
	doF(i, 40)
	{
		if(2*i < param.distThres) continue;
		cout<<"\r"<<i*2;
		param.distThres = _f 2 * i;
		VoteInitMask(); 
		cvS res = GetVoteError();
		fout2<<res.val[0]<<" "<<res.val[1]<<" "<<res.val[2]<<endl;
	}
	fout2.close();

}

string GPMAnalysisProc::GetShadowSegDir(string imgName)
{
	return "mask//" + imgName.substr(0, imgName.size()-4) + "s.png";
}
string GPMAnalysisProc::GetMaterialSegDir(string imgName)
{
	return "mask//" + imgName.substr(0, imgName.size()-4) + "m.png";
}


void GPMAnalysisProc::GetStatistics(ofstream& fout)
{
	GPMStatis sta;

	//cout<<"\nBegin analysis matching result for all images...\n";
	omp_set_num_threads(6);
#pragma omp parallel for
	doF(i, _i m_imgNames.size())
	{
		if(m_imgNames[i][m_imgNames[i].size()-5] == 'n') continue;
		//cout<<"--Analysis image "<<m_imgNames[i]<<"...\n";
		GPMStatis sta1;

		string corrFile = GetCorrFileDir(m_imgNames[i]);
		string shadowMaskFile = GetShadowSegDir(m_imgNames[i]);
		string matMaskDir = GetMaterialSegDir(m_imgNames[i]);

		cvi* shMaskOri = cvlig(m_fileDir + shadowMaskFile);
		cvi* shMask = cvci81(shMaskOri->width / param.downRatio, shMaskOri->height / param.downRatio);
		cvResize(shMaskOri, shMask, CV_INTER_NN);
		//cvsi(shMask);
		cvi* matMaskOri = cvlig(m_fileDir + matMaskDir);
		cvi* matMask = cvci(shMask);
		cvResize(matMaskOri, matMask, CV_INTER_NN);
		cvri(shMaskOri); cvri(matMaskOri);

		DenseCorrBox2D box;
		box.Load(m_fileDir + corrFile);

		doFcvi(shMask, i, j)
		{
			MultiCorr corrs = box.GetCorrsPerGrid(i, j);

			int shdV = _i cvg20(shMask, i, j), matV = _i cvg20(matMask, i, j);
			
			//unshadow
			if(shdV > 0)
			{
				doFv(k, corrs)
				{
					Corr& cor = corrs[k];
					sta1.unShdwMatchN ++;
					if(cor.dist <= param.distThres)
					{
						sta1.wum++;
					}
				}
				continue;
			}

			bool firstMatch = false;
			sta1.pixelN++;
			doFv(k, corrs)
			{
				Corr& cor = corrs[k];
				sta1.matchN++;
				if(cor.dist > param.distThres)
				{
					sta1.wm++;
					continue;
				}
				int shdV2 = _i cvg20(shMask, round(cor.x), round(cor.y));
				int matV2 = _i cvg20(matMask, round(cor.x), round(cor.y));
				if(shdV != shdV2 && matV == matV2)
				{
					if(firstMatch == false)
					{
						sta1.rp++;
						firstMatch = true;
					}
					sta1.dssm++;
				}
				else if(shdV != shdV2 && matV != matV2) sta1.dsdm++;
				else if(shdV == shdV2 && matV == matV2){
					//cout<<i<<" "<<j<<" "<<cor.x<<" "<<cor.y<<" "<<shdV<<" "<<shdV2<<" "<<cor.dist<<" "<<endl; pause;
					sta1.sssm++;
				}
				else sta1.ssdm++;
			}
		}
		//sta1.print();
		sta.add(sta1);

		cvri(shMask); cvri(matMask);
	}
	//cout<<"--Analysis matching result for all images complete! Summary: "<<endl;
	fout<<param.distThres<<"\t";
	sta.print2(fout);
	fout<<endl;
}

void GPMAnalysisProc::VoteInitMask()
{
	cout<<"Begin analysis voting result for all images...\n";

	omp_set_num_threads(6);
#pragma omp parallel for
	doF(ii, _i m_imgNames.size())
	{
		if(m_imgNames[ii][m_imgNames[ii].size()-5] == 'n') continue;
		cout<<"\rAnalysis image "<<m_imgNames[ii]<<"...";

		string corrFile = GetCorrFileDir(m_imgNames[ii]);
		string imgPrefix = m_imgNames[ii].substr(0, m_imgNames[ii].size() - 4);

		DenseCorrBox2D box;
		box.Load(m_fileDir + corrFile);

		ImgContainer img(m_fileDir + m_imgNames[ii], param.downRatio, param.colorMode);
		if(img.srcR() == NULL) img.GenerateResizedImg(1);

		vector<double> allGbV(img.src()->width * img.src()->height * 6);
		doFcvi(img.src(), i, j)
		{
			MultiCorr corrs = box.GetCorrsPerGrid(i, j);

			vector<double> gbV(6);

			Patch srcPatch, dstPatch;
			PatchDistMetric::GetPatch(img, _f i, _f j, 1.f, 0.f, (param.patchSize-1)/2, srcPatch);
			
			vector<vector<double> > gbVs(0);
			vector<pair<double, int> > Ls(0); 

			doFv(k, corrs)
			{
				Corr& v = corrs[k];
				if(v.dist > param.distThres) continue;
				PatchDistMetric::GetPatch(img, v.x, v.y, v.s, v.r, (param.patchSize-1)/2, dstPatch);
				double dist;
				double L = OptimizeGbVLab(srcPatch, dstPatch, gbV, dist);
				Ls.push_back(make_pair(L, Ls.size()));
				gbVs.push_back(gbV);
			}

			if(Ls.size() == 0)
			{
				doF(k, 3) gbV[2*k] = 1;
				doF(k, 3) gbV[2*k+1] = 0;
			}
			else
			{
				sort(Ls.begin(), Ls.end(), [](const pair<double, int>& v1, const pair<double, int>& v2){
					return v1.first > v2.first;
				});
				int idx = min2(1, _i Ls.size()-1);
				gbV = gbVs[Ls[idx].second];
			}

			doF(k, 6) allGbV[(i*img.src()->width+j)*6+k] = gbV[k];
		}

		doF(k, 3)
		{
			cvi* v1 = ShdwAnlysisProc::VislzVector(img.src()->width, img.src()->height, allGbV, 2*k, 0, 1);
			cvsi(m_fileDir + parVoteDir + imgPrefix + "_" + toStr(k) + "_gain.png", v1);
			cvi* v2 = ShdwAnlysisProc::VislzVector(img.src()->width, img.src()->height, allGbV, 2*k+1, -100, 100);
			cvsi(m_fileDir + parVoteDir + imgPrefix + "_" + toStr(k) + "_bias.png", v2);
			cvri(v1); cvri(v2);
			ShdwAnlysisProc::SaveVector(img.src()->width, img.src()->height, allGbV, 2*k, 
				m_fileDir + parVoteDir + imgPrefix + "_" + toStr(k) + "_gain.txt");
			ShdwAnlysisProc::SaveVector(img.src()->width, img.src()->height, allGbV, 2*k+1, 
				m_fileDir + parVoteDir + imgPrefix + "_" + toStr(k) + "_bais.txt");
		}

	}
	cout<<"\rVoting complete.\n";
}

cvS GPMAnalysisProc::GetVoteError()
{
	cout<<"Begin calculating voting error using GTParam and VoteParam...\n";
	cvS allError = cvs(0, 0, 0);
	doFv(i, m_imgNames)
	{
		if(m_imgNames[i][m_imgNames[i].size()-5] == 'n') continue;
		cout<<"\r--Analysis image "<<m_imgNames[i]<<"...";
		string imgPrefix = m_imgNames[i].substr(0, m_imgNames[i].size() - 4);

		string voteResDirLGain = m_fileDir + parVoteDir + imgPrefix + "_0_gain.txt";
		string gtResDirLGain = m_fileDir + parGTDir + imgPrefix + "_0_gain.txt";
		string voteResDirABias = m_fileDir + parVoteDir + imgPrefix + "_1_bais.txt";
		string gtResDirABias = m_fileDir + parGTDir + imgPrefix + "_1_bais.txt";
		string voteResDirBBias = m_fileDir + parVoteDir + imgPrefix + "_2_bais.txt";
		string gtResDirBBias = m_fileDir + parGTDir + imgPrefix + "_2_bais.txt";

		cvS sumError = cvs(0, 0, 0);
		ifstream fin(voteResDirLGain.c_str()), fin2(gtResDirLGain.c_str());
		ifstream fin3(voteResDirABias.c_str()), fin4(gtResDirABias.c_str());
		ifstream fin5(voteResDirBBias.c_str()), fin6(gtResDirBBias.c_str());
		int w, h, w2, h2;
		fin3>>w>>h; fin4>>w>>h; fin5>>w>>h; fin6>>w>>h; fin>>w>>h; fin2>>w2>>h2; 
		if(w != w2 || h != h2)
		{
			cout<<"Error!"<<endl;
		}
		else
		{
			double v1, v2, v3, v4, v5, v6;
			doF(k, w*h)
			{
				fin>>v1; fin2>>v2; fin3>>v3; fin4>>v4; fin5>>v5; fin6>>v6;
				sumError += cvs(abs(v1 - v2), abs(v3-v4), abs(v5-v6));
			}
			sumError /= w*h;
		}
		fin.close(); fin2.close();
		allError += sumError;
	}
	allError /= m_imgNames.size() / 2;
	cout<<"\rVote Error computing complete, vote error = "<<allError.val[0]<<" "<<allError.val[1]<<" "
		<<allError.val[2]<<".\n";
	return allError;
}

double GPMAnalysisProc::OptimizeGbVLab(Patch& srcPatch, Patch& dstPatch, vector<double>& gbV, double& dist)
{
	cvS srcAvg = cvs(0, 0, 0), dstAvg = cvs(0, 0, 0);
	doFv(k, srcPatch.hls)
	{
		srcAvg += srcPatch.hls[k];
	}
	srcAvg /= srcPatch.hls.size();
	doFv(k, dstPatch.hls)
	{
		dstAvg += dstPatch.hls[k];
	}
	dstAvg /= dstPatch.hls.size();
	
	gbV[0] = srcAvg.val[0] / dstAvg.val[0];
	gbV[1] = 0;
	gbV[2] = 1;
	gbV[3] = - srcAvg.val[1] + dstAvg.val[1];
	gbV[4] = 1;
	gbV[5] = - srcAvg.val[2] + dstAvg.val[2];

	cvS gain = cvs(gbV[0], gbV[2], gbV[4]);
	cvS bias = cvs(gbV[1], gbV[3], gbV[5]);
	float sum = 0;
	doFv(i, srcPatch.hls)
	{
		sum += _f cvSDSqr(srcPatch.hls[i], dstPatch.hls[i]*gain+bias);
	}
	dist = sqrt(sum / srcPatch.hls.size());

	return dstAvg.val[0];
}