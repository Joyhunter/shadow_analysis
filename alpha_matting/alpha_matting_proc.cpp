/*************************************************
  Copyright (C), 2014, Joy Artworks. Tsinghua Uni.
  All rights reserved.
  
  File name: alpha_matting_proc.cpp  v1.0
  Description: An  implementation of  "[07PAMI]  A 
               Closed-Form   Solution  to  Natural 
			   Image Matting"  , the most commonly
			   used alpha matting method.
  Other: Need OPENCV, sparse_mat_solver.

  Function List: 
    1. AlphaMattingProc::RunAlphaMatting
    2. s AlphaMattingProc::ParseErrorCode
    3. s AlphaMattingProc::Test

  History:
    1. Date: 2013.5.9
       Author: Li-Qian Ma
       Modification: Version 1.0.
*************************************************/
#include "stdafx.h"
#include "alpha_matting_proc.h"

namespace image_editing
{

	AlphaMattingProc::AlphaMattingProc(void){}
	AlphaMattingProc::~AlphaMattingProc(void){}

	int AlphaMattingProc::RunAlphaMatting(IN IplImage* src, IN IplImage* trimap, 
		OUT IplImage* back, OUT IplImage* fore, OUT IplImage* matte)
	{

			if(src == NULL || trimap == NULL || back == NULL || fore == NULL || matte == NULL)
			{
				return ALPHA_MATTING_NO_PRE_ALLOCATION;
			}

			if((cvGetSize(src) != cvGetSize(trimap)) || (cvGetSize(src) != cvGetSize(back))
				|| (cvGetSize(src) != cvGetSize(fore)) || (cvGetSize(src) != cvGetSize(matte)))
			{
				return ALPHA_MATTING_IMAGE_SIZE_NOT_EQUAL;
			}

			if((src->depth != IPL_DEPTH_8U) || (trimap->depth != IPL_DEPTH_8U) 
				|| (back->depth != IPL_DEPTH_8U) || (fore->depth != IPL_DEPTH_8U) 
				|| (matte->depth != IPL_DEPTH_8U))
			{
				return ALPHA_MATTING_IMAGE_DEPTH_NOT_MEET_REQUREMENT;
			}

			if((src->nChannels != 3) || (trimap->nChannels != 1) || (back->nChannels != 3)
				|| (fore->nChannels != 3) || (matte->nChannels != 1))
			{
				return ALPHA_MATTING_IMAGE_NCHANNEL_NOT_MEET_REQUREMENT;
			}

			bool has1 = false, has2 = false, has3 = false;
			doFcvi(trimap, i, j)
			{
				if(cvg20(trimap, i, j) == 255) has1 = true;
				else if(cvg20(trimap, i, j) == 0) has2 = true;
				else has3 = true;
			}
			if(!has1 || !has2 || !has3)
			{
				return ALPHA_MATTING_TRIMAP_NOT_MEET_REQUREMENT;
			}

			IplImage* tmpsrc, * tmpalp, * tmpfore, * tmpback, * tmptrimap;
			int x0,x1,y0,y1;
			x1 = y1 = 0; x0 = trimap->width-1; y0 = trimap->height-1;
			for(int i=0;i<trimap->height;i++){
				for(int j=0;j<trimap->width;j++){
					if (CV_IMAGE_ELEM(trimap,BYTE,i,j)){
						if (i>y1) y1 = i; if (i<y0) y0 = i;
						if (j>x1) x1 = j; if (j<x0) x0 = j;
					}
				}
			}
			//if(x1<x0||y1<y0) return;
			x0 = max(0,x0-5); y0 = max(0,y0-5); x1 = min(trimap->width-1,x1+5); y1 = min(trimap->height-1,y1+5);
			int w = x1-x0+1; int h = y1-y0+1;
			CvRect roi = cvRect(x0,y0,w,h);
			cvSetImageROI(src,roi);
			cvSetImageROI(trimap,roi);
			tmpsrc = cvCreateImage(cvSize(w,h),8,3);
			tmptrimap = cvCreateImage(cvSize(w,h),8,1);
			cvCopyImage(src,tmpsrc);
			cvCopyImage(trimap,tmptrimap);
			tmpalp = 0;
			tmpfore = 0;
			tmpback = 0;
			cvResetImageROI(src);
			cvResetImageROI(trimap);
			IplImage* alphaX;

			alphamatting::solveAlpha(tmpsrc, tmptrimap, tmpalp,alphaX);
			alphamatting::solveFBX(tmpsrc, alphaX, tmpfore, tmpback);

			cvReleaseImage(&alphaX);
			cvZero(fore);
			cvZero(matte);
			cvSetImageROI(matte,roi);
			cvSetImageROI(fore,roi);
			cvSetImageROI(back,roi);
			cvCopyImage(tmpalp,matte);
			cvCopyImage(tmpfore,fore);
			cvCopyImage(tmpback,back);
			cvResetImageROI(matte);
			cvResetImageROI(fore);
			cvResetImageROI(back);
			cvReleaseImage(&tmpalp);
			cvReleaseImage(&tmpback);
			cvReleaseImage(&tmpfore);
			cvReleaseImage(&tmpsrc);
			cvReleaseImage(&tmptrimap);
			
			return ALPHA_MATTING_SUCC;
	}

	string AlphaMattingProc::ParseErrorCode(int code)
	{
		switch(code)
		{
		case ALPHA_MATTING_SUCC:
			return "Alpha matting success.";
			break;
		case ALPHA_MATTING_NO_PRE_ALLOCATION:
			return "One or more of src/trimap/fore/back/alpha is not allocated memory before alpha matting.";
			break;
		case ALPHA_MATTING_IMAGE_SIZE_NOT_EQUAL:
			return "The size of src/trimap/fore/back/alpha is not totally equal.";
			break;
		case ALPHA_MATTING_IMAGE_DEPTH_NOT_MEET_REQUREMENT:
			return "The depth of src/trimap/fore/back/alpha doesn't totally meet requirement 8/8/8/8/8.";
			break;
		case ALPHA_MATTING_IMAGE_NCHANNEL_NOT_MEET_REQUREMENT:
			return "The nChannels of src/trimap/fore/back/alpha doesn't totally meet requirement 3/1/3/3/1.";
			break;
		case ALPHA_MATTING_TRIMAP_NOT_MEET_REQUREMENT:
			return "Trimap must have at least one pixel of foreground(255), background(0), and unknown(128).";
			break;
		default:
			return "Alpha matting invalid code.";
			break;
		}
	}

	void AlphaMattingProc::Test()
	{
		IplImage* src = cvLoadImage("alpha_matting_test_src.png", 1);
		IplImage* mask = cvLoadImage("alpha_matting_test_trimap.png", 0);
		if(!src || !mask){
			cout<<"'alpha_matting_test_src.png' or 'alpha_matting_test_trimap.png' is missing, exit!"<<endl;
			return;
		}

		IplImage* back = cvCloneImage(src);
		IplImage* fore = cvCloneImage(src);
		IplImage* matte = cvCloneImage(mask);

		cout<<"Alpha matting...\n";
		int errorCode = AlphaMattingProc().RunAlphaMatting(src, mask, back, fore, matte);
		cout<<ParseErrorCode(errorCode)<<endl;

		cvSaveImage("alpha_matting_test_fore.png", fore);
		cvSaveImage("alpha_matting_test_back.png", back);
		cvSaveImage("alpha_matting_test_matte.png", matte);

		cvReleaseImage(&src);
		cvReleaseImage(&mask);
		cvReleaseImage(&fore);
		cvReleaseImage(&back);
		cvReleaseImage(&matte);

	}

}