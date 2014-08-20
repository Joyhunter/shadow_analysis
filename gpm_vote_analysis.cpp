#include "StdAfx.h"
#include "gpm_analysis.h"
#include "shadow_analysis.h"

void GPMAnalysisProc::VoteAnalysis()
{
	wGetDirFiles(m_fileDir + "*.png", m_imgNames);

	param.downRatio = 6; param.multiLevelRatio = 2; param.extraLevelN = 0;
	param.patchSize = 7;
	param.colorMode = CV_BGR2Lab;

	rangeDir = "range3//";
	parVoteDir = "VoteParam//";
	parGTDir = "GTParam6//";

	TrainVoteRTrees();
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
	int matchN = 12, feaN = 4, yN = 3;
	int nPatches = 0;

	//compute nPatches
	doF(ii, _i m_imgNames.size())
	{
		if(m_imgNames[ii][m_imgNames[ii].size()-5] == 'n' || m_imgNames[ii][m_imgNames[ii].size()-5] == 'p') continue;
		cvi* src = cvlic(m_fileDir + m_imgNames[ii]);
		nPatches += src->width * src->height / param.downRatio / param.downRatio;
		cvri(src);
	}

	//init mat for rt
	Mat training_data = Mat(nPatches, matchN * feaN, CV_32FC1);  
	Mat training_target = Mat(nPatches, 1, CV_32FC1);

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
		cvi* gt = LoadGTFile(parGTDir, imgPrefix);

		vector<double> gbV(6);
		double dist;
		Patch srcPatch, dstPatch;
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

			//save feature to mat
			doFv(k, features)
			{
				training_data.at<float>(patchIdx, feaN*k) = features[k].lgain;
				training_data.at<float>(patchIdx, feaN*k+1) = features[k].abias;
				training_data.at<float>(patchIdx, feaN*k+2) = features[k].bbias;
				training_data.at<float>(patchIdx, feaN*k+3) = features[k].dist;
			}
			training_target.at<float>(patchIdx, 0) = _f cvg2(gt, i, j).val[0];
			patchIdx++;
		}

		cvri(gt);
	}
	cout<<"\nData ready.\n";

	//train
	Mat var_type = Mat(matchN * feaN + 1, 1, CV_8U);  
	var_type.setTo(Scalar(CV_VAR_NUMERICAL));
	//float priors[] = {1,1,1,1,1,1,1,1,1,1};  // weights of each classification for classes  
	CvRTParams params = CvRTParams(25,			// max depth  
									5,			// min sample count  
									0,			// regression accuracy: N/A here  
									false,		// compute surrogate split, no missing data  
									15,			// max number of categories (use sub-optimal algorithm for larger numbers)  
									NULL,		// the array of priors  
									false,		// calculate variable importance  
									4,			// number of variables randomly selected at node and used to find the best split(s).  
									100,		// max number of trees in the forest  
									0.01f,      // forrest accuracy  
									CV_TERMCRIT_ITER |   CV_TERMCRIT_EPS // termination cirteria  
									);
	cout<<"Training...";  
	CvRTrees* rtree = new CvRTrees;  
	rtree->train(training_data.rowRange(Range(0, 52500)), CV_ROW_SAMPLE, training_target.rowRange(Range(0, 52500)),  
		Mat(), Mat(), var_type, Mat(), params);
	cout<<"\rTraining complete!\n";
	CvFileStorage* fs = cvOpenFileStorage("example.xml", 0, CV_STORAGE_WRITE);
	rtree->write(fs, "value");

	//testing
	patchIdx = 0;
	float predictErrorSum = 0;
	doF(ii, _i m_imgNames.size())
	{
		if(m_imgNames[ii][m_imgNames[ii].size()-5] == 'n' || m_imgNames[ii][m_imgNames[ii].size()-5] == 'p') continue;
		string imgPrefix = m_imgNames[ii].substr(0, m_imgNames[ii].size() - 4);
		cvi* src = cvlic(m_fileDir + m_imgNames[ii]);
		cvi* res = cvci81(src->width / param.downRatio, src->height / param.downRatio);
		doFcvi(res, i, j)
		{
			Mat temp = training_data.row(patchIdx);
			float r = rtree->predict(temp);
			float target = training_target.at<float>(patchIdx, 0);
			predictErrorSum += fabs(r - target);
			cvs20(res, i, j, r*255);
			patchIdx++;
		}

		cvsi(m_fileDir + "TrainingRes//" + imgPrefix + "_predict.png", res);
		cvri(src); cvri(res);
	}
	cout<<"Predict error in training set: "<<predictErrorSum<<endl;

}

cvi* GPMAnalysisProc::LoadGTFile(string gtDir, string imgPrefix)
{
	string fileDir[3];
	fileDir[0] = m_fileDir + gtDir + imgPrefix + "_0_gain.txt";
	fileDir[1] = m_fileDir + gtDir + imgPrefix + "_1_bias.txt";
	fileDir[2] = m_fileDir + gtDir + imgPrefix + "_2_bias.txt";

	int w, h;
	ifstream fin1(fileDir[0].c_str());
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