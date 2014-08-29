#include "StdAfx.h"
#include "gpm_analysis.h"
#include "shadow_analysis.h"

void GPMAnalysisProc::InitRTParam()
{
	ifstream fin("rt.cfg");
	fin>>param.downRatio>>param.multiLevelRatio>>param.extraLevelN>>param.patchSize;
	fin>>rangeDir;
	fin>>parGTDir;

	fin>>rtParam.matchNPerPatch;
	fin>>rtParam.sampleStep;
	fin>>rtParam.max_depth;
	fin>>rtParam.min_sample_count;
	fin>>rtParam.reggression_accurancy;
	fin>>rtParam.max_categories;
	fin>>rtParam.nactive_vars;
	fin>>rtParam.max_tree_N;
	fin>>rtParam.forest_accurancy;
	fin>>rtParam.forest_accurancyAB;
	fin>>parPredictDir;
	fin.close();
}

void GPMAnalysisProc::VoteAnalysis()
{
	wGetDirFiles(m_fileDir + "*.png", m_imgNames);

	param.downRatio = 6; param.multiLevelRatio = 2; param.extraLevelN = 0;
	param.patchSize = 7;
	param.colorMode = CV_BGR2Lab;

	rangeDir = "range4//";
	//parVoteDir = "VoteParam//";
	parGTDir = "GTParam6//";
	parPredictDir = "PredictParam//";

	rtParam.matchNPerPatch = 12;
	rtParam.sampleStep = 1;
	rtParam.max_depth = 25;
	rtParam.min_sample_count = 5;
	rtParam.reggression_accurancy = 0.f;
	rtParam.max_categories = 15;
	rtParam.nactive_vars = 4;
	rtParam.max_tree_N = 100;
	rtParam.forest_accurancy = 0.01f;
	rtParam.forest_accurancyAB = 1.f;
		
	InitRTParam();

	TrainVoteRTrees();

	//PredictUseRTrees();
}

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
	int matchN = rtParam.matchNPerPatch, feaN = 4, yN = 3;
	int nPatches = 0;

	//compute nPatches
	doF(ii, _i m_imgNames.size())
	{
		if(m_imgNames[ii][m_imgNames[ii].size()-5] == 'n' || m_imgNames[ii][m_imgNames[ii].size()-5] == 'p') continue;
		cvi* src = cvlic(m_fileDir + m_imgNames[ii]);
		nPatches += _i ceil(_f(src->width / param.downRatio) / rtParam.sampleStep) * 
			_i ceil(_f(src->width / param.downRatio) / rtParam.sampleStep);
		cvri(src);
	}
	cout<<"Npatch = "<<nPatches<<". FeatureSize = "<<matchN * feaN<<".\n";

	//init mat for rt
	Mat training_data = Mat(nPatches, matchN * feaN, CV_32FC1);  
	Mat training_target = Mat(nPatches, 1, CV_32FC1);
	Mat training_target2 = Mat(nPatches, 1, CV_32FC1);
	Mat training_target3 = Mat(nPatches, 1, CV_32FC1);

	//get data to mat
	int patchIdx = 0;
	cout<<"Get GPM result for all images...\n";
	doF(ii, _i m_imgNames.size())
	{
		if(m_imgNames[ii][m_imgNames[ii].size()-5] == 'n' || m_imgNames[ii][m_imgNames[ii].size()-5] == 'p') continue;
		cout<<"\rAnalysis image "<<m_imgNames[ii]<<"...";

		//load image and coor file
		string corrFile = GetCorrFileDir(m_imgNames[ii]);
		string imgPrefix = m_imgNames[ii].substr(0, m_imgNames[ii].size() - 4);
		DenseCorrBox2D box;
		box.Load(m_fileDir + corrFile);
		ImgContainer img(m_fileDir + m_imgNames[ii], param.downRatio, param.colorMode);
		if(img.srcR() == NULL) img.GenerateResizedImg(1);
		
		//load groundtruth output
		cvi* gt = LoadGTFile(m_fileDir + parGTDir, imgPrefix);

		vector<double> gbV(6);
		double dist;
		Patch srcPatch, dstPatch;
		doFcvi(img.src(), i, j)
		{
			if(i % rtParam.sampleStep != 0 || j % rtParam.sampleStep != 0) continue;

			PatchDistMetric::GetPatch(img, _f i, _f j, 1.f, 0.f, (param.patchSize-1)/2, srcPatch);

			MultiCorr corrs = box.GetCorrsPerGrid(i, j);
			vector<MatchFeature> features(corrs.size());
			doFv(k, corrs)
			{
				Corr& v = corrs[k];
				features[k].LoadFromCorr(v);
				//recompute gain and bias
				PatchDistMetric::GetPatch(img, v.x, v.y, v.s, v.r, (param.patchSize-1)/2, dstPatch);
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
	cout<<"\nData ready.\n";

	//train
	Mat var_type = Mat(matchN * feaN + 1, 1, CV_8U);  
	var_type.setTo(Scalar(CV_VAR_NUMERICAL));
	//float priors[] = {1,1,1,1,1,1,1,1,1,1};  // weights of each classification for classes  
	CvRTParams params = CvRTParams(rtParam.max_depth,			// max depth  
									rtParam.min_sample_count,			// min sample count  
									rtParam.reggression_accurancy,			// regression accuracy: N/A here  
									false,		// compute surrogate split, no missing data  
									rtParam.max_categories,			// max number of categories (use sub-optimal algorithm for larger numbers)  
									NULL,		// the array of priors  
									false,		// calculate variable importance  
									rtParam.nactive_vars,			// number of variables randomly selected at node and used to find the best split(s).  
									rtParam.max_tree_N,		// max number of trees in the forest  
									rtParam.forest_accurancy,      // forrest accuracy  
									CV_TERMCRIT_ITER |   CV_TERMCRIT_EPS // termination cirteria  
									);
	cout<<"Training...";  
	CvFileStorage* fs = cvOpenFileStorage((m_fileDir + "rt.xml").c_str(), 0, CV_STORAGE_WRITE);

	CvRTrees* rtree = new CvRTrees;  
	rtree->train(training_data, CV_ROW_SAMPLE, training_target,  
		Mat(), Mat(), var_type, Mat(), params);
	rtree->write(fs, "LGAIN");
	delete rtree;
	rtree = new CvRTrees;  
	CvRTParams params2 = CvRTParams(rtParam.max_depth, rtParam.min_sample_count, rtParam.reggression_accurancy,	false,
		rtParam.max_categories,	NULL, false, rtParam.nactive_vars, rtParam.max_tree_N, rtParam.forest_accurancyAB, 
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

	cout<<"\rTraining complete!\n";
}

void GPMAnalysisProc::PredictUseRTrees()
{
	cout<<"Read random forest file...\n";
	CvFileStorage* fs = cvOpenFileStorage((m_fileDir + "rt.xml").c_str(), 0, CV_STORAGE_READ);
	CvRTrees* rtree = new CvRTrees, *rtree2 = new CvRTrees, *rtree3 = new CvRTrees;
	rtree->read(fs, cvGetFileNodeByName(fs, 0, "LGAIN")); cout<<1;
	rtree2->read(fs, cvGetFileNodeByName(fs, 0, "ABIAS"));
	rtree3->read(fs, cvGetFileNodeByName(fs, 0, "BBIAS"));
	cvReleaseFileStorage(&fs);

	wMkDir(m_fileDir + parPredictDir);

	int matchN = rtParam.matchNPerPatch, feaN = 4;
	Mat dataRow = Mat(1, matchN * feaN, CV_32FC1);

	cout<<"Fetch GPM result...\n";
	cvS sumAll = cvs(0, 0, 0);
	doF(ii, _i m_imgNames.size())
	{
		if(m_imgNames[ii][m_imgNames[ii].size()-5] == 'n' || m_imgNames[ii][m_imgNames[ii].size()-5] == 'p') continue;
		cout<<"Analysis image "<<m_imgNames[ii]<<"...";

		//load image and coor file
		string corrFile = GetCorrFileDir(m_imgNames[ii]);
		string imgPrefix = m_imgNames[ii].substr(0, m_imgNames[ii].size() - 4);
		DenseCorrBox2D box;
		box.Load(m_fileDir + corrFile);
		ImgContainer img(m_fileDir + m_imgNames[ii], param.downRatio, param.colorMode);
		if(img.srcR() == NULL) img.GenerateResizedImg(1);

		//load groundtruth output
		cvi* gt = LoadGTFile(m_fileDir + parGTDir, imgPrefix);

		vector<double> gbV(6);
		double dist;
		Patch srcPatch, dstPatch;
		vector<double> predictV(6*img.src()->width*img.src()->height);
		cvS errorSum = cvs(0, 0, 0);
		doFcvi(img.src(), i, j)
		{
			PatchDistMetric::GetPatch(img, _f i, _f j, 1.f, 0.f, (param.patchSize-1)/2, srcPatch);

			MultiCorr corrs = box.GetCorrsPerGrid(i, j);
			vector<MatchFeature> features(corrs.size());
			doFv(k, corrs)
			{
				Corr& v = corrs[k];
				features[k].LoadFromCorr(v);
				//recompute gain and bias
				PatchDistMetric::GetPatch(img, v.x, v.y, v.s, v.r, (param.patchSize-1)/2, dstPatch);
				double L = OptimizeGbVLab(srcPatch, dstPatch, gbV, dist);
				features[k].lgain = _f gbV[0];
				features[k].abias = _f gbV[3];
				features[k].bbias = _f gbV[5];
				features[k].dist = _f dist;
			}
			sort(features.begin(), features.end(), [](const MatchFeature& v1, const MatchFeature& v2){
				return (v1.lgain == v2.lgain)?(v1.dist < v2.dist):(v1.lgain < v2.lgain);
			});

			doF(k, matchN) //doFv(k, features)
			{
				dataRow.at<float>(0, feaN*k) = features[k].lgain;
				dataRow.at<float>(0, feaN*k+1) = features[k].abias;
				dataRow.at<float>(0, feaN*k+2) = features[k].bbias;
				dataRow.at<float>(0, feaN*k+3) = features[k].dist;
			}
			float r1 = rtree->predict(dataRow), r2 = rtree2->predict(dataRow), r3 = rtree3->predict(dataRow);
			predictV[(i*img.src()->width+j)*6] = r1;
			predictV[(i*img.src()->width+j)*6+3] = r2;
			predictV[(i*img.src()->width+j)*6+5] = r3;

			//error
			cvS v = cvg2(gt, i, j);
			errorSum.val[0] += abs(v.val[0] - _d r1);
			errorSum.val[1] += abs(v.val[1] - _d r2);
			errorSum.val[2] += abs(v.val[2] - _d r3);
		}
		errorSum /= img.src()->width*img.src()->height;
		cout<<"\nError = "<<errorSum.val[0]<<","<<errorSum.val[1]<<","<<errorSum.val[2]<<"\n";
		sumAll += errorSum;

		//save predictV
		doF(k, 3)
		{
			cvi* v1 = ShdwAnlysisProc::VislzVector(img.src()->width, img.src()->height, predictV, 2*k, 0, 1);
			cvsi(m_fileDir + parPredictDir + imgPrefix + "_" + toStr(k) + "_gain.png", v1);
			cvi* v2 = ShdwAnlysisProc::VislzVector(img.src()->width, img.src()->height, predictV, 2*k+1, -100, 100);
			cvsi(m_fileDir + parPredictDir + imgPrefix + "_" + toStr(k) + "_bias.png", v2);
			cvri(v1); cvri(v2);
			ShdwAnlysisProc::SaveVector(img.src()->width, img.src()->height, predictV, 2*k, 
				m_fileDir + parPredictDir + imgPrefix + "_" + toStr(k) + "_gain.txt");
			ShdwAnlysisProc::SaveVector(img.src()->width, img.src()->height, predictV, 2*k+1, 
				m_fileDir + parPredictDir + imgPrefix + "_" + toStr(k) + "_bais.txt");
		}

		cvri(gt);
	}
	cout<<"\nComplete.\n";

	sumAll /= m_imgNames.size() / 2;
	cout<<"\nError = "<<sumAll.val[0]<<","<<sumAll.val[1]<<","<<sumAll.val[2]; pause;
}

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
}