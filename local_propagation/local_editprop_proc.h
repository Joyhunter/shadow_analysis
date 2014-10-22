/*************************************************
  Copyright (C), 2014, Joy Artworks. Tsinghua Uni.
  All rights reserved.
  
  File name: local_editprop_proc.h  v1.0
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

/*************************************************
  MACRO: SOLVE_SPARSE_SYSTEM_USING_MATLAB
  Description: solve this sparse linear system 
               using matlab or C++ implementation
    - Matlab   need matlab_engine.h/lib
	           add matlab include/lib dir to proj.
	- C++      need sparse_mat_solver header
	           less robust than Matlab.
*************************************************/
//#define SOLVE_SPARSE_SYSTEM_USING_MATLAB

#ifdef SOLVE_SPARSE_SYSTEM_USING_MATLAB
#include "../matlab_engine/matlab_engine.h"
#pragma comment(lib, "../Release/matlab_engine.lib")
#else
#include "../matrix/matrix_engine.h"
#endif

namespace image_editing
{

#define LOCAL_EDITPROP_SUCC 0
#define LOCAL_EDITPROP_NO_PRE_ALLOCATION 1
#define LOCAL_EDITPROP_IMAGE_SIZE_NOT_EQUAL 2
#define LOCAL_EDITPROP_IMAGE_DEPTH_NOT_MEET_REQUREMENT 3
#define LOCAL_EDITPROP_IMAGE_NCHANNEL_NOT_MEET_REQUREMENT 4
#define LOCAL_EDITPROP_TRIMAP_NOT_MEET_REQUREMENT 5

	class LocalEditpropProc
	{

	public:

		LocalEditpropProc(void);
		~LocalEditpropProc(void);

		/*************************************************
		  Function: RunManifoldEditprop
		  Description: Run [06TOG] local edit prop.
		  Calls: None
		  Called By: ::Test()
		  Input: src:    source image. allocated 8U3C.
		         stroke: stroke image. allocated 8U1C.
				         0   - user constrain not modify
						 128 - others
						 255 - user interest region
				 alpha:  referred in paper. 0.1f ~ 5.0f
				         larger makes prop easier.
				 eps:    referred in paper. epsilon
				 nambda: referred in paper. 0.0002f ~ 20.0
				         larger makes prop easier.
		  Output: result: result edit mask. allocated 8U1C
		  Return: Error code.
		  Others: result must be allocated before calling
		          size(src) == size(stroke) ==size(result)
				  depth(src/stroke/result) = 8
				  nChannels(src/stroke/result) = 3/1/1
		*************************************************/
		int RunLocalEditprop(IN cvi* src, IN cvi* stroke, 
			OUT cvi* result, 
			IN float alpha = 0.5f, IN float eps = 0.0001f, 
			IN float nambda = 1.0f);		

		/*************************************************
		  Function: ParseErrorCode
		  Description: Parse return information of 
		               RunLocalEditprop.
		  Calls: None
		  Called By: ::Test().
		  Input: code:    return value of RunLocalEditprop.
		  Output: String that explains the return value.
		*************************************************/
		static string ParseErrorCode(int code);

		/*************************************************
		  Function: Test
		  Description: Run local edit propagation test. 
		  Calls: ::RunLocalEditprop()
		  Called By: None.
		  Other: 'editprop_test_src.png' and 
		         'editprop_test_stroke.png'  is read 
				 from disk. Resulting edit mask
				 'local_editprop_test_result.png', is
				 saved to disk.
		*************************************************/
		static void Test();

	private:



	};

}