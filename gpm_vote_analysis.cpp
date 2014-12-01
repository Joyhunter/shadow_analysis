#include "StdAfx.h"
#include "gpm_analysis.h"
#include "shadow_analysis.h"
#include "param_loader.h"

struct MatchFeature
{
	float lgain, abias, bbias, dist;
	float x, y;
	void LoadFromCorr(Corr& corr)
	{
		lgain = _f corr.gain.val[0];
		abias = _f corr.bias.val[1];
		bbias = _f corr.bias.val[2];
		dist = _f corr.dist;
		x = corr.x; y = corr.y;
	}
};

//(1)2211.79 (3)2116.86 (4)2122.1 (6)4091.96 (9)2122.32
void GPMAnalysisProc::TrainVoteRTrees()
{
	int matchN = vcfg.rtParam.matchNPerPatch, feaN = 4, yN = 3;
	int nPatches = 0;

	COutYel("-- Random Forest Training --\n");
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

	//compute nPatches
	COutGreen("Step 1: Calculate all patches...\n");
	doFv(idx, m_imgNames)
	{
		string fileName = m_imgNames[idx];
		string imgPrefix = m_imgNames[idx].substr(0, m_imgNames[idx].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		if(vcfg.focusedPrefix.size() > 0 && find(vcfg.focusedPrefix.begin(), vcfg.focusedPrefix.end(), imgPrefix) == vcfg.focusedPrefix.end())
			continue;
		COutTeal("\rProcessing " + imgPrefix + ".png..");

		string corrFile = vcfg.srcDir + vcfg.corrDir + imgPrefix + ".corr";
		DenseCorrBox2D box;
		box.Load(corrFile);
		if(box.m_dstH == 0)
		{
			COutRed("Coor file " + corrFile + " not exist error.\n");
			continue;
		}
		int w = box.m_dstW, h = box.m_dstH;

		nPatches += _i ceil(_f w / vcfg.rtParam.sampleStep) * _i ceil(_f h / vcfg.rtParam.sampleStep);
	}
	COutGreen("\rNpatch = " + toStr(nPatches) + ". FeatureSize = " + toStr(matchN * feaN) + ".\n");

	//init mat for rt
	Mat training_data = Mat(nPatches, matchN * feaN, CV_32FC1);  
	Mat training_target = Mat(nPatches, 1, CV_32FC1);
	Mat training_target2 = Mat(nPatches, 1, CV_32FC1);
	Mat training_target3 = Mat(nPatches, 1, CV_32FC1);

	//get data to mat
	int patchIdx = 0;
	COutGreen("Step 2: Get Data Ready...\n");
	doFv(idx, m_imgNames)
	{
		string fileName = m_imgNames[idx];
		string imgPrefix = m_imgNames[idx].substr(0, m_imgNames[idx].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		if(vcfg.focusedPrefix.size() > 0 && find(vcfg.focusedPrefix.begin(), vcfg.focusedPrefix.end(), imgPrefix) == vcfg.focusedPrefix.end())
			continue;
		COutTeal("\rProcessing " + imgPrefix + ".png..");

		//load image and coor file
		string corrFile = vcfg.srcDir + vcfg.corrDir + imgPrefix + ".corr";
		DenseCorrBox2D box;
		box.Load(corrFile);
		if(box.m_dstH == 0)
		{
			COutRed("Coor file " + corrFile + " not exist error.\n");
			continue;
		}
		int w = box.m_dstW, h = box.m_dstH;
		ImgContainer img(vcfg.srcDir + fileName, 1.0, cfg.src_colorMode_gpm);
		img.GenerateResizedImg(w, h);
		
		//load groundtruth output
		cvi* gt = ParamLoader::LoadParamFromDir(vcfg.srcDir + vcfg.gtSubDir, imgPrefix);

		vector<double> gbV(6);
		double dist;
		Patch srcPatch, dstPatch;
		doFcvi(img.srcR(), i, j)
		{
			if(i % vcfg.rtParam.sampleStep != 0 || j % vcfg.rtParam.sampleStep != 0) continue;

			PatchDistMetric::GetPatch(img, _f i, _f j, 1.f, 0.f, false, false, (vcfg.patchSize-1)/2, srcPatch);

			MultiCorr corrs = box.GetCorrsPerGrid(i, j);
			vector<MatchFeature> features(corrs.size());
			doFv(k, corrs)
			{
				Corr& v = corrs[k];
				features[k].LoadFromCorr(v);
				//recompute gain and bias
				PatchDistMetric::GetPatch(img, v.x, v.y, v.s, v.r, v.hr, v.vr, (vcfg.patchSize-1)/2, dstPatch);
				double L = OptimizeGbVLab(srcPatch, dstPatch, gbV, dist);
				features[k].lgain = _f gbV[0];
				features[k].abias = _f gbV[3];
				features[k].bbias = _f gbV[5];
				features[k].dist = _f dist;
			}
			sort(features.begin(), features.end(), [](const MatchFeature& v1, const MatchFeature& v2){
				return (v1.lgain == v2.lgain)?(v1.dist < v2.dist):(v1.lgain < v2.lgain);
			});

			//save feature to mat
			while(_i features.size() < matchN) features.push_back(features[features.size()-1]);
			doF(k, matchN) //doFv(k, features)
			{
				training_data.at<float>(patchIdx, feaN*k) = features[k].lgain;
				training_data.at<float>(patchIdx, feaN*k+1) = features[k].abias;
				training_data.at<float>(patchIdx, feaN*k+2) = features[k].bbias;
				training_data.at<float>(patchIdx, feaN*k+3) = features[k].dist;
			}
			training_target.at<float>(patchIdx, 0) = _f cvg2(gt, i, j).val[0];
			training_target2.at<float>(patchIdx, 0) = _f cvg2(gt, i, j).val[1];
			training_target3.at<float>(patchIdx, 0) = _f cvg2(gt, i, j).val[2];
			patchIdx++;
		}

		cvri(gt);
	}
	COutGreen("\rStep 2: Data Ready...\n");

	//train
	COutGreen("Step 3: Training...\n");
	Mat var_type = Mat(matchN * feaN + 1, 1, CV_8U);  
	var_type.setTo(Scalar(CV_VAR_NUMERICAL));
	//float priors[] = {1,1,1,1,1,1,1,1,1,1};  // weights of each classification for classes  
	CvRTParams params = CvRTParams(vcfg.rtParam.max_depth,			// max depth  
									vcfg.rtParam.min_sample_count,			// min sample count  
									vcfg.rtParam.reggression_accurancy,			// regression accuracy: N/A here  
									false,		// compute surrogate split, no missing data  
									vcfg.rtParam.max_categories,			// max number of categories (use sub-optimal algorithm for larger numbers)  
									NULL,		// the array of priors  
									false,		// calculate variable importance  
									vcfg.rtParam.nactive_vars,			// number of variables randomly selected at node and used to find the best split(s).  
									vcfg.rtParam.max_tree_N,		// max number of trees in the forest  
									vcfg.rtParam.forest_accurancy,      // forrest accuracy  
									CV_TERMCRIT_ITER |   CV_TERMCRIT_EPS // termination cirteria  
									);
	CvFileStorage* fs = cvOpenFileStorage((vcfg.srcDir + vcfg.modelFile).c_str(), 0, CV_STORAGE_WRITE);

	CvRTrees* rtree = new CvRTrees;  
	rtree->train(training_data, CV_ROW_SAMPLE, training_target,  
		Mat(), Mat(), var_type, Mat(), params);
	rtree->write(fs, "LGAIN");
	delete rtree;
	rtree = new CvRTrees;  
	CvRTParams params2 = CvRTParams(vcfg.rtParam.max_depth, vcfg.rtParam.min_sample_count, vcfg.rtParam.reggression_accurancy,	false,
		vcfg.rtParam.max_categories,	NULL, false, vcfg.rtParam.nactive_vars, vcfg.rtParam.max_tree_N, vcfg.rtParam.forest_accurancyAB, 
		CV_TERMCRIT_ITER |   CV_TERMCRIT_EPS);
	rtree->train(training_data, CV_ROW_SAMPLE, training_target2,  
		Mat(), Mat(), var_type, Mat(), params2);
	rtree->write(fs, "ABIAS");
	delete rtree;
	rtree = new CvRTrees;  
	rtree->train(training_data, CV_ROW_SAMPLE, training_target3,  
		Mat(), Mat(), var_type, Mat(), params2);
	rtree->write(fs, "BBIAS");
	delete rtree;

	cvReleaseFileStorage(&fs);

	COutGreen("Step 3: Training done.\n");
}


void GPMAnalysisProc::PredictUseRTrees()
{
	COutYel("-- Random Forest Prediction --\n");
	wMkDir(vcfg.srcDir + vcfg.resSubdir);
	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);

	//Read random forest file
	COutGreen("Reading rt model file..\n");
	CvFileStorage* fs = cvOpenFileStorage((vcfg.srcDir + vcfg.modelFile).c_str(), 0, CV_STORAGE_READ);
	CvRTrees* rtree = new CvRTrees, *rtree2 = new CvRTrees, *rtree3 = new CvRTrees;
	rtree->read(fs, cvGetFileNodeByName(fs, 0, "LGAIN"));
	rtree2->read(fs, cvGetFileNodeByName(fs, 0, "ABIAS"));
	rtree3->read(fs, cvGetFileNodeByName(fs, 0, "BBIAS"));
	cvReleaseFileStorage(&fs);

	int matchN = vcfg.rtParam.matchNPerPatch, feaN = 4;
	Mat dataRow = Mat(1, matchN * feaN, CV_32FC1);

	cvS sumAll = cvs(0, 0, 0);
	doFv(idx, m_imgNames)
	{
		string fileName = m_imgNames[idx];
		string imgPrefix = m_imgNames[idx].substr(0, m_imgNames[idx].size() - 4);
		if(imgPrefix.size() >= 2 && imgPrefix.substr(imgPrefix.size()-2, 2) == "_n") continue;
		if(vcfg.focusedPrefix.size() > 0 && find(vcfg.focusedPrefix.begin(), vcfg.focusedPrefix.end(), imgPrefix) == vcfg.focusedPrefix.end())
			continue;
		COutTeal("Processing " + imgPrefix + ".png..\n");

		//load image and coor file
		string corrFile = vcfg.srcDir + vcfg.corrDir + imgPrefix + ".corr";
		DenseCorrBox2D box;
		box.Load(corrFile);
		if(box.m_dstH == 0)
		{
			COutRed("Coor file " + corrFile + " not exist error.\n");
			continue;
		}
		int w = box.m_dstW, h = box.m_dstH;
		ImgContainer img(vcfg.srcDir + fileName, 1.0, cfg.src_colorMode_gpm);
		img.GenerateResizedImg(w, h);

		//load groundtruth output
		//cvi* gt = LoadGTFile(cfg.srcDir + parGTDir, imgPrefix);

		vector<double> gbV(6);
		double dist;
		Patch srcPatch, dstPatch;
		cvi* param = cvci323(w, h);
		cvS errorSum = cvs(0, 0, 0);
		doFcvi(img.srcR(), i, j)
		{
			PatchDistMetric::GetPatch(img, _f i, _f j, 1.f, 0.f, false, false, (vcfg.patchSize-1)/2, srcPatch);

			MultiCorr corrs = box.GetCorrsPerGrid(i, j);
			vector<MatchFeature> features(corrs.size());
			doFv(k, corrs)
			{
				Corr& v = corrs[k];
				features[k].LoadFromCorr(v);
				//recompute gain and bias
				PatchDistMetric::GetPatch(img, v.x, v.y, v.s, v.r, v.hr, v.vr, (vcfg.patchSize-1)/2, dstPatch);
				double L = OptimizeGbVLab(srcPatch, dstPatch, gbV, dist);
				features[k].lgain = _f gbV[0];
				features[k].abias = _f gbV[3];
				features[k].bbias = _f gbV[5];
				features[k].dist = _f dist;
			}
			sort(features.begin(), features.end(), [](const MatchFeature& v1, const MatchFeature& v2){
				return (v1.lgain == v2.lgain)?(v1.dist < v2.dist):(v1.lgain < v2.lgain);
			});

			while(_i features.size() < matchN) features.push_back(features[features.size()-1]);
			doF(k, matchN)
			{
				dataRow.at<float>(0, feaN*k) = features[k].lgain;
				dataRow.at<float>(0, feaN*k+1) = features[k].abias;
				dataRow.at<float>(0, feaN*k+2) = features[k].bbias;
				dataRow.at<float>(0, feaN*k+3) = features[k].dist;
			}

			float r1 = rtree->predict(dataRow), r2 = rtree2->predict(dataRow), r3 = rtree3->predict(dataRow);
			cvs2(param, i, j, cvs(r1, r2, r3));

			//error
			//cvS v = cvg2(gt, i, j);
			//errorSum.val[0] += abs(v.val[0] - _d r1);
			//errorSum.val[1] += abs(v.val[1] - _d r2);
			//errorSum.val[2] += abs(v.val[2] - _d r3);
		}
		//errorSum /= img.src()->width*img.src()->height;
		//cout<<"\nError = "<<errorSum.val[0]<<","<<errorSum.val[1]<<","<<errorSum.val[2]<<"\n";
		//sumAll += errorSum;

		//save predictV
		ParamLoader::SaveParamToDir(vcfg.srcDir + vcfg.resSubdir, imgPrefix, param);
		ParamLoader::ShowParamInDir(vcfg.srcDir + vcfg.resSubdir, imgPrefix, param);
		cvri(param);

		//cvri(gt);
	}

	//sumAll /= m_imgNames.size() / 2;
	//cout<<"\nError = "<<sumAll.val[0]<<","<<sumAll.val[1]<<","<<sumAll.val[2]; pause;
}

/*
cvi* GPMAnalysisProc::LoadGTFile(string gtDir, string imgPrefix)
{
	string fileDir[3];
	fileDir[0] = gtDir + imgPrefix + "_0_gain.txt";
	fileDir[1] = gtDir + imgPrefix + "_1_bais.txt";
	fileDir[2] = gtDir + imgPrefix + "_2_bais.txt";

	int w, h;
	ifstream fin1(fileDir[0].c_str());
	if(!fin1)
	{
		cout<<"LoadGTFile: "<<fileDir[0]<<" not exist! error!"<<endl;
		return 0;
	}
	fin1>>w>>h;
	cvi* res = cvci323(w, h);

	ifstream fin2(fileDir[1].c_str());
	ifstream fin3(fileDir[2].c_str());
	fin2>>w>>h; fin3>>w>>h;
	double temp[3];
	doFcvi(res, i, j)
	{
		fin1>>temp[0]; fin2>>temp[1]; fin3>>temp[2];
		cvs2(res, i, j, cvs(temp[0], temp[1], temp[2]));
	}

	fin1.close(); fin2.close(); fin3.close();
	return res;
	}*/

void GPMAnalysisProc::InitRTParam()
{
	// 	ifstream fin("rt.cfg");
	// 	fin>>param.downRatio>>param.multiLevelRatio>>param.extraLevelN>>param.patchSize;
	// 	fin>>rangeDir;
	// 	fin>>parGTDir;
	// 
	// 	fin>>rtParam.matchNPerPatch;
	// 	fin>>rtParam.sampleStep;
	// 	fin>>rtParam.max_depth;
	// 	fin>>rtParam.min_sample_count;
	// 	fin>>rtParam.reggression_accurancy;
	// 	fin>>rtParam.max_categories;
	// 	fin>>rtParam.nactive_vars;
	// 	fin>>rtParam.max_tree_N;
	// 	fin>>rtParam.forest_accurancy;
	// 	fin>>rtParam.forest_accurancyAB;
	// 	fin>>parPredictDir;
	// 	fin.close();
}

void GPMAnalysisProc::VoteAnalysis()
{
	// 	wGetDirFiles(cfg.srcDir + "*.png", m_imgNames);
	// 
	// 	param.downRatio = 6; param.multiLevelRatio = 2; param.extraLevelN = 0;
	// 	param.patchSize = 7;
	// 	param.colorMode = CV_BGR2Lab;
	// 
	// 	rangeDir = "range4//";
	// 	//parVoteDir = "VoteParam//";
	// 	parGTDir = "GTParam6//";
	// 	parPredictDir = "PredictParam//";
	// 
	// 	rtParam.matchNPerPatch = 12;
	// 	rtParam.sampleStep = 1;
	// 	rtParam.max_depth = 25;
	// 	rtParam.min_sample_count = 5;
	// 	rtParam.reggression_accurancy = 0.f;
	// 	rtParam.max_categories = 15;
	// 	rtParam.nactive_vars = 4;
	// 	rtParam.max_tree_N = 100;
	// 	rtParam.forest_accurancy = 0.01f;
	// 	rtParam.forest_accurancyAB = 1.f;
	// 		
	// 	InitRTParam();
	// 
	// 	TrainVoteRTrees();
	// 
	// 	//PredictUseRTrees();
}
