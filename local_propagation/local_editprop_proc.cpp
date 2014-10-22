/*************************************************
  Copyright (C), 2014, Joy Artworks. Tsinghua Uni.
  All rights reserved.
  
  File name: local_editprop_proc.cpp  v1.0
  Description: An direct implementation of  "
               [06TOG]Interactive Local Adjustment
			   of Tonal Value". Options(C++/matlab 
			   sparse solver) can be modified.
  Other: Need OPENCV, sparse_mat(C++ solver)
         / matlab_eng(Matlab solver).

  Macro: 
    1. SOLVE_SPARSE_SYSTEM_USING_MATLAB
	   use C++ solver - light (default)
	   use Matlab solver - robust

  Function List: 
    1. LocalEditpropProc::RunLocalEditprop
	2. s LocalEditpropProc::ParseErrorCode
	3. s LocalEditpropProc::Test

  History:
    1. Date: 2014.5.23
       Author: Li-Qian Ma
       Modification: Version 1.0.
*************************************************/
#pragma once
#include "stdafx.h"
#include "local_editprop_proc.h"

namespace image_editing
{

	LocalEditpropProc::LocalEditpropProc(void){}

	LocalEditpropProc::~LocalEditpropProc(void){}

	void GetEditsFromStroke(cvi* stroke, vecI& vEditIndices, vecF& vEditValues)
	{
		vEditIndices.clear(); vEditValues.clear();
		doFcvi(stroke, i, j)
		{
			float v = (float)cvg20(stroke, i, j);
			if(v == 128) continue;
			vEditIndices.push_back(i * stroke->width + j);
			vEditValues.push_back(v / 255.0f);
		}
	}

	void SolveLinearSystem(IN cvi* src, IN vecI& vEditIndices, IN vecF& vEditValues,
		OUT vecD& propRes,
		float alpha, float eps, float nambda)
	{
		int nrows = src->width*src->height;
		vecI I, J;
		vecD F;
		vecD b(nrows);

		int editIdx = 0, idx = 0;
		doFcvi(src, i, j){

			double w = 0, g = 0;
			if(vEditIndices.size() > (size_t)editIdx && idx == vEditIndices[editIdx]){
				w = 1;
				g = vEditValues[editIdx];
				editIdx ++;
			}

			double sum = 0;
			for (int oi = -1; oi <= 1; oi ++) for (int oj = -1; oj <= 1; oj ++){
				if(abs(oi) + abs(oj) != 1) continue;
				if(cvIn(i+oi, j+oj, src)){
					double A = cvSD(cvg2(src, i, j), cvg2(src, i+oi, j+oj)) / 255.0f;
					A = - nambda / (pow((float)A, alpha) + eps);
					I.push_back(idx); J.push_back((i+oi)*src->width+j+oj);
					F.push_back(A);
					sum += A;
				}
			}

			I.push_back(idx); J.push_back(idx);
			F.push_back(w - sum);
			b[idx] = w * g;

			idx ++;
		}

		propRes.resize(nrows);
#ifdef SOLVE_SPARSE_SYSTEM_USING_MATLAB
		Utility::MatFunction().
			solveSparseLinearSystem(nrows, nrows, I, J, F, b, propRes);
#else
		Utility::MatrixFunction().
			solveSparseLinearSystem(nrows, nrows, I, J, F, b, propRes);
#endif

	}

	int LocalEditpropProc::RunLocalEditprop(IN cvi* src, IN cvi* stroke, 
		OUT cvi* result, 
		IN float alpha, IN float eps, IN float nambda)
	{

		if(src == NULL || stroke == NULL || result == NULL)
			return LOCAL_EDITPROP_NO_PRE_ALLOCATION;

		if((cvGetSize(src) != cvGetSize(stroke)) || (cvGetSize(src) != cvGetSize(result)))
			return LOCAL_EDITPROP_IMAGE_SIZE_NOT_EQUAL;

		if((src->depth != IPL_DEPTH_8U) || (stroke->depth != IPL_DEPTH_8U) 
			|| (result->depth != IPL_DEPTH_8U))
			return LOCAL_EDITPROP_IMAGE_DEPTH_NOT_MEET_REQUREMENT;

		if((src->nChannels != 3) || (stroke->nChannels != 1) || (result->nChannels != 1))
			return LOCAL_EDITPROP_IMAGE_NCHANNEL_NOT_MEET_REQUREMENT;

		bool has1 = false, has2 = false;
		doFcvi(stroke, i, j){
			if(cvg20(stroke, i, j) > 128) has1 = true;
			else if(cvg20(stroke, i, j) == 0) has2 = true;
		}
		if(!has1 || !has2)
			return LOCAL_EDITPROP_TRIMAP_NOT_MEET_REQUREMENT;

		vecI vEditIndices; vecF vEditValues;
		GetEditsFromStroke(stroke, vEditIndices, vEditValues);

		cvi* srcLab = cvci(src);
		cvCvtColor(src, srcLab, CV_BGR2Lab);

		vecD propRes;
		SolveLinearSystem(srcLab, vEditIndices, vEditValues, propRes, alpha, eps, nambda);
		
		doFcvi(result, i, j)
		{
			cvs20(result, i, j, propRes[i*result->width + j] * 255.0f);
		}

		return LOCAL_EDITPROP_SUCC;
	}

	string LocalEditpropProc::ParseErrorCode(int code)
	{
		switch(code)
		{
		case LOCAL_EDITPROP_SUCC:
			return "Local edit propagation success.";
			break;
		case LOCAL_EDITPROP_NO_PRE_ALLOCATION:
			return "One or more of src/stroke/result is not allocated memory before alpha matting.";
			break;
		case LOCAL_EDITPROP_IMAGE_SIZE_NOT_EQUAL:
			return "The size of src/stroke/result is not totally equal.";
			break;
		case LOCAL_EDITPROP_IMAGE_DEPTH_NOT_MEET_REQUREMENT:
			return "The depth of src/stroke/result doesn't totally meet requirement 8/8/8.";
			break;
		case LOCAL_EDITPROP_IMAGE_NCHANNEL_NOT_MEET_REQUREMENT:
			return "The nChannels of src/stroke/result doesn't totally meet requirement 3/1/1.";
			break;
		case LOCAL_EDITPROP_TRIMAP_NOT_MEET_REQUREMENT:
			return "Trimap must have at least one pixel of user interest(255), non modified(0).";
			break;
		default:
			return "Local edit propagation invalid code.";
			break;
		}
	}

	void LocalEditpropProc::Test()
	{
		cvi* src = cvlic("editprop_test_src.png");
		cvi* stroke = cvlig("editprop_test_stroke.png");
		if(!src || !stroke){
			cout<<"'editprop_test_src.png' or 'editprop_test_stroke.png' is missing, exit!"<<endl;
			return;
		}

		cout<<"Local edit propagation...\n";
		cvi* result = cvci(stroke);
		int errorCode = LocalEditpropProc().RunLocalEditprop(src, stroke, result);
		cout<<ParseErrorCode(errorCode)<<endl;

		cvsi("local_editprop_test_result.png", result);
		cvri(src);
		cvri(stroke);
		cvri(result);
	}
}