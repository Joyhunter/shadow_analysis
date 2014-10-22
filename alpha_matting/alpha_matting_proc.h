/*************************************************
  Copyright (C), 2014, Joy Artworks. Tsinghua Uni.
  All rights reserved.
  
  File name: alpha_matting_proc.h  v1.0
  Description: An  implementation of  "[07PAMI]  A 
               Closed-Form   Solution  to  Natural 
			   Image Matting" ,  the most commonly
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
#pragma once

#include "include/global.h"

#ifdef _DEBUG
#pragma comment(lib, "sparsematd.lib")
#else
#pragma comment(lib, "sparsemat.lib")
#endif

namespace image_editing
{

#define ALPHA_MATTING_SUCC 0
#define ALPHA_MATTING_NO_PRE_ALLOCATION 1
#define ALPHA_MATTING_IMAGE_SIZE_NOT_EQUAL 2
#define ALPHA_MATTING_IMAGE_DEPTH_NOT_MEET_REQUREMENT 3
#define ALPHA_MATTING_IMAGE_NCHANNEL_NOT_MEET_REQUREMENT 4
#define ALPHA_MATTING_TRIMAP_NOT_MEET_REQUREMENT 5

	class AlphaMattingProc
	{

	public:

		AlphaMattingProc(void);
		~AlphaMattingProc(void);

		/*************************************************
		  Function: RunAlphaMatting
		  Description: Run [07PAMI] alpha matting.
		  Calls: ./include/*.cpp
		  Called By: ::Test()
		  Input: img:    source image. allocated 8U3C.
		         triMap: trimap image. allocated 8U1C.
				         0   - background
						 128 - unknown
				         255 - foreground 
		  Output: back:	 allocated 8U3C. 
		                 background colors for each pixel.
				  fore:	 allocated 8U3C. 
						 foreground colors for each pixel.
				  alpha: allocated 8U1C. 
				         alpha channel. 
		  Return: Error code.
		  Others: img = fore * alpha + back * (1 - alpha).
		          img, trimap, back, fore,  alpha  must be
				  allocated before the call of this  func.
				  fore & alpha can be directly  import  to 
				  photoshop to do further design.
		*************************************************/
		int RunAlphaMatting(IN IplImage* img, IN IplImage* triMap, 
			IN OUT IplImage* back, IN OUT IplImage* fore, IN OUT IplImage* alpha);

		/*************************************************
		  Function: ParseErrorCode
		  Description: Parse return information of 
		               RunAlphaMatting.
		  Calls: None
		  Called By: ::Test().
		  Input: code:    return value of RunAlphaMatting.
		  Output: String that explains the return value.
		*************************************************/
		static string ParseErrorCode(int code);

		/*************************************************
		  Function: Test
		  Description: Run alpha matting test. 
		  Calls: ::RunAlphaMatting()
		  Called By: None.
		  Other: 'alpha_matting_test_src.png' and 
		         'alpha_matting_test_trimap.png' is read 
				 from disk. Alpha matting results
				 'alpha_matting_test_fore.png',
				 'alpha_matting_test_back.png',
				 'alpha_matting_test_alpha.png' is saved
				 to disk.
		*************************************************/
		static void Test();
	};

}