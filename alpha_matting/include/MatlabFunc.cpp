#include "stdafx.h"
#include "global.h"

#define epsilon		(1e-7)

namespace alphamatting {

	void GetMean(CvMat *src, CvMat *dst) {
		CvSize size = cvGetSize(src);
		for (int i = 0; i < size.width; i++) {
			double mean = 0;
			for (int j = 0; j < size.height; j++)
				mean += cvmGet(src, j, i);
			cvmSet(dst, i, 0, mean / size.height);
		}
	}

	void repmat(CvMat *src, int height, int width, CvMat *dst) {
		CvSize size = cvGetSize(src);
		for (int r = 0; r < size.height * height; r++)
			for (int c = 0; c < size.width * width; c++)
				cvmSet(dst, r, c, cvmGet(src, r % size.height, c % size.width));
	}

	int qsortcompare(const void* elem1, const void* elem2)  {
		int e1 = *(int *)elem1;
		int e2 = *(int *)elem2;
		if (e1 < e2)
			return -1;
		else if (e1 == e2)
			return 0;
		else if (e1 > e2)
			return 1;
		return 0;
	}

	/*
	输入：
		img			原图
		consts		已知值的区域。
						0表示为原trimap中的未知区域
						1表示为原trimap中的已知区域
		win_size	未知区域膨胀的窗口半径大小
	输出：
		A			Laplace矩阵	
	*/
	void GetLaplacian(IplImage *img, IplImage *consts, int win_size, CvSparseMat *&A, int *&nonZeroIndex, int *&allNonZero) {
		if (A)
			cvReleaseSparseMat(&A);

		//neb_size=(win_size*2+1)^2;
		int neb_size = (win_size * 2 + 1) * (win_size * 2 + 1);

		//[h,w,c]=size(I);
		int h = img->height;
		int w = img->width;
		int c = img->nChannels;

		//n=h; m=w;
		int n = h;
		int m = w;

		//img_size=w*h;
		int img_size = w * h;

		//consts=imerode(consts,ones(win_size*2+1));	

		IplImage *tmpconsts = cvCloneImage(consts);
		////for (int cr = 0; cr < consts->height; cr++)
		////	for (int cc = 0; cc < consts->width; cc++)
		////		if (CV_IMAGE_ELEM(consts, BYTE, cr, cc) == 0) {
		////			RECT rect;
		////			GetRect(consts, cr, cc, win_size, &rect);
		////			for (int rr = rect.top; rr <= rect.bottom; rr++)
		////				for (int rc = rect.left; rc <= rect.right; rc++)						
		////					CV_IMAGE_ELEM(tmpconsts, BYTE, rr, rc) = 0;
		////		}
		cvErode(tmpconsts, tmpconsts);
		consts = tmpconsts;
		//cvErode(consts, consts);				//未知区域膨胀一圈

		////PM01(consts, "consts.txt");
		////exit(0);
		////ok..print(consts, "consts.txt");	

		//indsM=reshape([1:img_size],h,w);
		CvMat *indsM = cvCreateMat(h, w, CV_64FC1);	
		for (int ic = 0; ic < w; ic++)
			for (int ir = 0; ir < h; ir++)
				cvmSet(indsM, ir, ic, ic * h + ir);

		//tlen=sum(sum(1-consts(win_size+1:end-win_size,win_size+1:end-win_size)))*(neb_size^2);
		int tlen = 0;
		for (int tr = 1; tr < consts->height - 1; tr++)
			for (int tc = 1; tc < consts->width - 1; tc++)
				tlen += 1 - CV_IMAGE_ELEM(consts, BYTE, tr, tc);		//	计算去掉边框1个宽度未知区域的大小
		tlen *= neb_size * neb_size;									//	计算总的未知数个数				

		//row_inds=zeros(tlen ,1);
		CvMat *row_inds = cvCreateMat(tlen, 1, CV_64FC1);
		cvZero(row_inds);				//	未知数所在的行下标
		
		//col_inds=zeros(tlen,1);
		CvMat *col_inds = cvCreateMat(tlen, 1, CV_64FC1);
		cvZero(col_inds);				//	未知数所在的列下标

		//vals=zeros(tlen,1);
		CvMat *vals = cvCreateMat(tlen, 1, CV_64FC1);
		cvZero(vals);					//	未知数的系数

		//len=0;
		int len = 0;

		CvMat *win_inds = cvCreateMat(neb_size, 1, CV_64FC1);
		CvMat *winI = cvCreateMat(neb_size, c, CV_64FC1);
		CvMat *win_mu = cvCreateMat(c, 1, CV_64FC1);
		CvMat *win_var = cvCreateMat(c, c, CV_64FC1);
		CvMat *tvals = cvCreateMat(neb_size, neb_size, CV_64FC1);

		CvMat *tmp1 = cvCreateMat(c, c, CV_64FC1);
		CvMat *tmp2 = cvCreateMat(c, c, CV_64FC1);
		CvMat *tmp3 = cvCreateMat(neb_size, c, CV_64FC1);
		CvMat *tmp4 = cvCreateMat(c, neb_size, CV_64FC1);
		CvMat *tmp5 = cvCreateMat(neb_size * 1, 1 * neb_size, CV_64FC1);
		CvMat *tmp6 = cvCreateMat(1 * neb_size, neb_size * 1, CV_64FC1);
		CvMat *tmp7 = cvCreateMat(1, neb_size, CV_64FC1);

		//for j=1+win_size:w-win_size
		for (int j = win_size; j < w - win_size; j++) {
		//	for i=win_size+1:h-win_size
			for (int i = win_size; i < h - win_size; i++) {
		//		if (consts(i,j))
				if (CV_IMAGE_ELEM(consts, BYTE, i, j)) {
		//			continue
					continue;
		//		end  
				}

				cvZero(win_inds);
				cvZero(winI);
				cvZero(win_mu);
				cvZero(win_var);
				cvZero(tvals);

				cvZero(tmp1);
				cvZero(tmp2);
				cvZero(tmp3);
				cvZero(tmp4);
				cvZero(tmp5);
				cvZero(tmp6);
				cvZero(tmp7);

				RECT rect;
				rect.left = j - win_size;
				rect.right = j + win_size;
				rect.top = i - win_size;
				rect.bottom = i + win_size;

		//		win_inds=indsM(i-win_size:i+win_size,j-win_size:j+win_size);
		//		win_inds=win_inds(:);						
				GetRect(indsM, &rect, win_inds);		//	计算以未知像素为中心的窗口中每一个像素的下标
				//PM(win_inds);

		//		winI=I(i-win_size:i+win_size,j-win_size:j+win_size,:);
		//		winI=reshape(winI,neb_size,c);			
				GetRect(img, &rect, c, neb_size, winI);	//	获得原图中的值
														//	winI: 9x3
				////printf("winI = \n");
				////PM(winI);

		//		win_mu=mean(winI,1)';			
				GetMean(winI, win_mu);					//	计算未知像素的每一个通道上的平均值
														//	win_mu:	1x3
				////printf("win_mu = \n");
				////PM(win_mu);
				
		//		win_var=inv(winI'*winI/neb_size-win_mu*win_mu' +epsilon/neb_size*eye(c));						
				cvMulTransposed(winI, tmp1, 1);			//	tmp1: 3x3
				cvMulTransposed(win_mu, tmp2, 0);		//	tmp2: 3x3

				////printf("winI' * winI = \n");
				////PM(tmp1);

				////printf("win_mu * win_mu' = \n");
				////PM(tmp2);
				for (int varr = 0; varr < c; varr++)
					for (int varc = 0; varc < c; varc++)  {
						double value = cvmGet(tmp1, varr, varc) / neb_size - cvmGet(tmp2, varr, varc) + epsilon / neb_size * (varr == varc);
						cvmSet(win_var, varr, varc, value);
					}
				
				////PM(win_var);
				cvInvert(win_var, win_var);
				////PM(win_var);
				////exit(0);

		//		winI=winI-repmat(win_mu',neb_size,1);
				for (int winr = 0; winr < neb_size; winr++)
					for (int winc = 0; winc < c; winc++) {
						double value = cvmGet(winI, winr, winc) - cvmGet(win_mu, winc, 0);
						cvmSet(winI, winr, winc, value);
					}
				//PM(winI);

		//		tvals=(1+winI*win_var*winI')/neb_size;			
				cvGEMM(winI, win_var, 1, NULL, 0, tmp3);
				cvTranspose(winI, tmp4);
				cvGEMM(tmp3, tmp4, 1, NULL, 0, tvals);
				for (int r = 0; r < neb_size; r++)
					for (int c = 0; c < neb_size; c++) {
						double v = cvmGet(tvals, r, c);
						cvmSet(tvals, r, c, (v + 1) / neb_size);
					}
				//PM(tvals, "tvals.txt");

		//		row_inds(1+len:neb_size^2+len)=reshape(repmat(win_inds,1,neb_size),...
		//			neb_size^2,1);		
				repmat(win_inds, 1, neb_size, tmp5);
				for (int p = 0; p < neb_size * neb_size; p++) {
					double v = cvmGet(tmp5, p % neb_size, p / neb_size);
					cvmSet(row_inds, p + len, 0, v);				
				}	

		//		col_inds(1+len:neb_size^2+len)=reshape(repmat(win_inds',neb_size,1),...
		//			neb_size^2,1);
				cvTranspose(win_inds, tmp7);
				repmat(tmp7, neb_size, 1, tmp6);
				for (int p = 0; p < neb_size * neb_size; p++) {
					double v = cvmGet(tmp6, p % neb_size, p / neb_size);
					cvmSet(col_inds, p + len, 0, v);
				}

		//		vals(1+len:neb_size^2+len)=tvals(:);
				for (int p = 0; p < neb_size * neb_size; p++) {
					double v = cvmGet(tvals, p % neb_size, p / neb_size);
					cvmSet(vals, p + len, 0, v);
				}

		//		len=len+neb_size^2;
				len += neb_size * neb_size;
		//	end
			}
		//end  
		} 
		//SHOWTIME("Coefficient");

		//vals=vals(1:len);
		//row_inds=row_inds(1:len);
		//col_inds=col_inds(1:len);
		//A=sparse(row_inds,col_inds,vals,img_size,img_size);
		//sumA=sum(A,2);
		//A=spdiags(sumA(:),0,img_size,img_size)-A;

		nonZeroIndex = new int[img_size];
		memset(nonZeroIndex, 255, _msize(nonZeroIndex));
		allNonZero = new int[tlen];
		memset(allNonZero, 0, _msize(allNonZero));

		CvMat *sumA = cvCreateMat(img_size, 1, CV_64FC1);
		cvZero(sumA);
		for (int i = 0; i < len; i++) {
			int tr = (int)(cvmGet(row_inds, i, 0) + 0.5);
			nonZeroIndex[tr] = 0;
			double v = cvmGet(sumA, tr, 0);
			v += cvmGet(vals, i, 0);
			cvmSet(sumA, tr, 0, v);
		}
		//SHOWTIME("Get Sum A");

		int nonZeroTotal = 0;
		for (int i = 0; i < img_size; i++) {
			if (nonZeroIndex[i] == 0) {
				nonZeroIndex[i] = nonZeroTotal;
				allNonZero[nonZeroTotal] = i;
				nonZeroTotal++;			
			}
		}
		//SHOWTIME("Build Nonzero Index");

		int Asize[2] = {nonZeroTotal, nonZeroTotal};
		A = cvCreateSparseMat(2, Asize, CV_64FC1);
		//SHOWTIME("Alloc Mat");

		for (int i = 0; i < len; i++) {
			int tr = (int)(cvmGet(row_inds, i, 0) + 0.5);
			int tc = (int)(cvmGet(col_inds, i, 0) + 0.5);
			tr = nonZeroIndex[tr];
			tc = nonZeroIndex[tc];
			double v = cvGetReal2D(A, tr, tc);
			v -= cvmGet(vals, i, 0);
			cvSetReal2D(A, tr, tc, v);
		}

		for (int i = 0; i < nonZeroTotal; i++) {
			double v = cvmGet(sumA, allNonZero[i], 0) + cvGetReal2D(A, i, i);
			cvSetReal2D(A, i, i, v);
		}
		//SHOWTIME("Get -A");		

		////for (int i = 0; i < len; i++) {
		////	int tr = cvmGet(row_inds, i, 0) + 0.5;
		////	int tc = cvmGet(col_inds, i, 0) + 0.5;
		////	double v = cvGetReal2D(A, tr, tc);
		////	v -= cvmGet(vals, i, 0);
		////	cvSetReal2D(A, tr, tc, v);
		////}
		////SHOWTIME("Get -A");

		//////int *allrow = new int[len];
		//////for (int i = 0; i < len; i++) {
		//////	allrow[i] = cvmGet(row_inds, i, 0) + 0.5;
		//////}
		//////qsort(allrow, len, sizeof(int), qsortcompare);
		//////for (int i = 0; i < len; i++) {
		//////	if (i == 0 || (i != 0 && allrow[i] != allrow[i - 1])) {
		//////		double v = cvmGet(sumA, allrow[i], 0) + cvGetReal2D(A, allrow[i], allrow[i]);
		//////		cvSetReal2D(A, allrow[i], allrow[i], v);
		//////	}
		//////}
		//////delete[] allrow;

		//PM(A, "qsortA.txt");

		////for (int i = 0; i < img_size; i++) {
		////	double v = cvmGet(sumA, i, 0) + cvGetReal2D(A, i, i);
		////	cvSetReal2D(A, i, i, v);
		////}
		////PM(A, "dirA.txt");
		//////SHOWTIME("Set Mat Value");
		cvReleaseMat(&win_inds);
		cvReleaseMat(&winI);
		cvReleaseMat(&win_mu);
		cvReleaseMat(&win_var);
		cvReleaseMat(&tvals);

		cvReleaseMat(&tmp1);
		cvReleaseMat(&tmp2);
		cvReleaseMat(&tmp3);
		cvReleaseMat(&tmp4);
		cvReleaseMat(&tmp5);
		cvReleaseMat(&tmp6);
		cvReleaseMat(&tmp7);

		cvReleaseMat(&indsM);
		cvReleaseMat(&row_inds);
		cvReleaseMat(&col_inds);
		cvReleaseMat(&vals);	
		cvReleaseMat(&sumA);

		cvReleaseImage(&tmpconsts);
		//SHOWTIME("Release All");
	}

	//function alpha=solveAlpha(I,consts_map,consts_vals,varargin)
	void solveAlpha(IplImage *src, IplImage *trimap, IplImage *&alpha, IplImage *&alpha01) {
		IplImage *consts_map = cvCreateImage(cvGetSize(trimap), 8, 1);
		cvZero(consts_map);
		
		for (int r = 0; r < trimap->height; r++)
			for (int c = 0; c < trimap->width; c++) {
				if (CV_IMAGE_ELEM(trimap, BYTE, r, c) > 0.9 * 255 ||
					CV_IMAGE_ELEM(trimap, BYTE, r, c) < 0.1 * 255)
					CV_IMAGE_ELEM(consts_map, BYTE, r, c) = 1;
			}
		
		IplImage *consts_vals = cvCreateImage(cvGetSize(trimap), 8, 1);
		cvZero(consts_vals);
		for (int r = 0; r < trimap->height; r++)
			for (int c = 0; c < trimap->width; c++) {
				if (CV_IMAGE_ELEM(trimap, BYTE, r, c) > 0.9 * 255)				
					CV_IMAGE_ELEM(consts_vals, BYTE, r, c) = 1;
			}

		//A=getLaplacian1(I,consts_map,varargin{:});
		CvSparseMat *A = NULL;
		int *nonZeroIndex = NULL;
		int *allNonZero = NULL;		
		//SHOWTIME("Build consts");
		GetLaplacian(src, consts_map, 1, A, nonZeroIndex, allNonZero);
		
		//SHOWTIME("Get Laplacian");
		//print(A, "input\\sing\\BeforeA.txt");
		
		//[h,w,c]=size(I);	
		int h = src->height;
		int w = src->width;
		int c = src->nChannels;

		//img_size=w*h;	
		int img_size = w * h;

		//D=spdiags(consts_map(:),0,img_size,img_size);	
		//lambda=100;	
		int lambda = 100;
		//x=(A+lambda*D)\(lambda*consts_map(:).*consts_vals(:));

		int totalNonZero = cvGetDimSize(A, 0);
		for (int i = 0; i < totalNonZero; i++) {
			double v = cvGetReal2D(A, i, i);
			int tr = allNonZero[i] % h;
			int tc = allNonZero[i] / h;
			cvSetReal2D(A, i, i, v + lambda * CV_IMAGE_ELEM(consts_map, BYTE, tr, tc));
		}
		
		//print(A, "input\\sing\\NewA.txt");
					
		CvMat *b = cvCreateMat(totalNonZero, 1, CV_64FC1);
		CvMat *x = cvCreateMat(totalNonZero, 1, CV_64FC1);		
		for (int i = 0; i < totalNonZero; i++) {
			int tr = allNonZero[i] % h;
			int tc = allNonZero[i] / h;
			double v = lambda * CV_IMAGE_ELEM(consts_map, BYTE, tr, tc) * CV_IMAGE_ELEM(consts_vals, BYTE, tr, tc);
			cvmSet(b, i, 0, v);
		}

		//printf("*******************before\n");
		//SHOWTIME();
		//SolveSymSparseMat(A, b, x, 2, 1e-5);

		KunXuSolver(A, b, x);
		
		//sparse_matrix::SolveSymSparseMatUnsymmetric(A, b, x);
		//printf("*******************after\n");
		//SHOWTIME();
		//print(x, "input\\sing\\NewX.txt");
		
		//alpha=max(min(reshape(x,h,w),1),0);	
		if (alpha)
			cvReleaseImage(&alpha);

		alpha01 = cvCreateImage(cvGetSize(trimap), IPL_DEPTH_64F, 1);
		for (int r = 0; r < trimap->height; r++) {
			for (int c = 0; c < trimap->width; c++) {
				if (CV_IMAGE_ELEM(trimap, BYTE, r, c) == 255) {
					CV_IMAGE_ELEM(alpha01, double, r, c) = 1.0;
				}
				else if (CV_IMAGE_ELEM(trimap, BYTE, r, c) == 0) {
					CV_IMAGE_ELEM(alpha01, double, r, c) = 0.0;
				}
			}
		}

		alpha = cvCloneImage(trimap);
		for (int i = 0; i < totalNonZero; i++) {
			int tr = allNonZero[i] % h;
			int tc = allNonZero[i] / h;
			double v = cvmGet(x, i, 0);
			v = MAX(MIN(v, 1), 0);
			CV_IMAGE_ELEM(alpha, BYTE, tr, tc) = (int)(floor(v * 255 + 0.5));
			CV_IMAGE_ELEM(alpha01, double, tr, tc) = v;
		}

		cvReleaseSparseMat(&A);
		cvReleaseMat(&b);
		cvReleaseMat(&x);
		delete[] nonZeroIndex;
		delete[] allNonZero;
	}

	//function I = imIndexToVect(Y,X,imHeight)
	int imIndexToVector(int Y, int X, int imHeight) {
		//I = reshape(Y + (X-1)*imHeight,prod(size(X)),1); 
		return Y + X * imHeight;
	}

	void sparse(CvMat *inds1, CvMat *inds2, CvMat *vals, int total, int height, int width, CvSparseMat *&m) {
		if (m)
			cvReleaseSparseMat(&m);
		int sizes[2] = {height, width};
		m = cvCreateSparseMat(2, sizes, CV_64FC1);
		for (int i = 0; i < total; i++) {
			cvSetReal2D(m, (int)cvmGet(inds1, i, 0), (int)cvmGet(inds2, i, 0), cvmGet(vals, i, 0));
		}
	}

	//function [Gx,Gy,G3,G4]=getGMatByMask(w,h,mask)
	void GetGMatByMask(int w, int h, IplImage *mask, CvSparseMat *&Gx, CvSparseMat *&Gy, CvSparseMat *&G3, CvSparseMat *&G4) {
		//imgSize=w*h;
		int imgSize = w * h;

		//dS=[1,-1];
		CvMat *dS = cvCreateMat(1, 2, CV_64FC1);
		cvmSet(dS, 0, 0, 1);
		cvmSet(dS, 0, 1, -1);

		//filtSizeS=1;
		int filtSizeS = 1;

		//%indsGx1=[]; indsGx2=[]; valsGx=[];
		//%indsGy1=[]; indsGy2=[]; valsGy=[];
		//indsGx1=zeros(imgSize*2,1);
		CvMat *indsGx1 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsGx1);

		//indsGx2=zeros(imgSize*2,1);
		CvMat *indsGx2 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsGx2);

		//valsGx=zeros(imgSize*2,1);
		CvMat *valsGx = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(valsGx);

		//indsGy1=zeros(imgSize*2,1);
		CvMat *indsGy1 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsGy1);

		//indsGy2=zeros(imgSize*2,1);
		CvMat *indsGy2 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsGy2);

		//valsGy=zeros(imgSize*2,1);
		CvMat *valsGy = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(valsGy);

		//indsG31=zeros(imgSize*2,1);
		CvMat *indsG31 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsG31);

		//indsG32=zeros(imgSize*2,1);
		CvMat *indsG32 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsG32);

		//valsG3=zeros(imgSize*2,1);
		CvMat *valsG3 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(valsG3);

		//indsG41=zeros(imgSize*2,1);
		CvMat *indsG41 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsG41);

		//indsG42=zeros(imgSize*2,1);
		CvMat *indsG42 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(indsG42);

		//valsG4=zeros(imgSize*2,1);
		CvMat *valsG4 = cvCreateMat(imgSize * 2, 1, CV_64FC1);
		cvZero(valsG4);

		//indy=0; indx=0; ind3=0; ind4=0;
		int indy = -1, indx = -1, ind3 = -1, ind4 = -1;

		//for x=1:w-1,	
		for (int x = 0; x < w - 1; x++) {
		//	for y=1:h,
			for (int y = 0; y < h; y++) {
		//		if ((~mask(y,x))&(~mask(y,x+1)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y, x + 1))) {
		//			continue
					continue;
		//		end  
				}
		//		for disp=0:filtSizeS,
				for (int disp = 0; disp <= filtSizeS; disp++) {
		//			indx=indx+1;
					indx = indx + 1;
		//			indsGx1(indx)=imIndexToVect(y,x,h);				
					cvmSet(indsGx1, indx, 0, imIndexToVector(y, x, h));
		//			indsGx2(indx)=imIndexToVect(y,x+disp,h);
					cvmSet(indsGx2, indx, 0, imIndexToVector(y, x + disp, h));
		//			valsGx(indx)=dS(disp+1);
					cvmSet(valsGx, indx, 0, cvmGet(dS, 0, disp));
		//		end
				}
		//	end
			}
		//end
		}
		
		//for x=1:w,	
		for (int x = 0; x < w; x++) {		
		//	for y=1:h-1,
			for (int y = 0; y < h - 1; y++) {
		//		if ((~mask(y,x))&(~mask(y+1,x)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y + 1, x))) {
		//			continue
					continue;
		//		end
				}
		//		for disp=0:filtSizeS,
				for (int disp = 0; disp <= filtSizeS; disp++) {
		//			indy=indy+1;
					indy = indy + 1;
		//			indsGy1(indy)=imIndexToVect(y,x,h);
					cvmSet(indsGy1, indy, 0, imIndexToVector(y, x, h));
		//			indsGy2(indy)=imIndexToVect(y+disp,x,h);
					cvmSet(indsGy2, indy, 0, imIndexToVector(y + disp, x, h));
		//			valsGy(indy)=dS(disp+1);
					cvmSet(valsGy, indy, 0, cvmGet(dS, 0, disp));
		//		end;
				}
		//	end;
			}
		//end
		}

		//for x=1:w-1,	
		for (int x = 0; x < w - 1; x++) {		
		//	for y=1:h-1,
			for (int y = 0; y < h - 1; y++) {
		//		if ((~mask(y,x))&(~mask(y+1,x+1)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y + 1, x + 1))) {
		//			continue
					continue;
		//		end
				}
		//		for disp=0:filtSizeS,
				for (int disp = 0; disp <= filtSizeS; disp++) {
		//			ind3=ind3+1;
					ind3 = ind3 + 1;
		//			indsG31(ind3)=imIndexToVect(y,x,h);
					cvmSet(indsG31, ind3, 0, imIndexToVector(y, x, h));
		//			indsG32(indy)=imIndexToVect(y+disp,x+disp,h);
					cvmSet(indsG32, ind3, 0, imIndexToVector(y + disp, x + disp, h));
		//			valsG3(ind3)=dS(disp+1);
					cvmSet(valsG3, ind3, 0, cvmGet(dS, 0, disp));
		//		end
				}
		//	end
			}			
		//end
		}

		//for x=1:w-1,	
		for (int x = 0; x < w - 1; x++) {
		//	for y=2:h,
			for (int y = 1; y < h; y++) {
		//		if ((~mask(y,x))&(~mask(y-1,x+1)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y - 1, x + 1))) {
		//			continue
					continue;
		//		end
				}
		//		for disp=0:filtSizeS,   
				for (int disp = 0; disp <= filtSizeS; disp++) {
		//			ind4=ind4+1;
					ind4 = ind4 + 1;
		//			indsG41(ind4)=imIndexToVect(y,x,h);
					cvmSet(indsG41, ind4, 0, imIndexToVector(y, x, h));
		//			indsG42(ind4)=imIndexToVect(y-disp,x+disp,h);
					cvmSet(indsG42, ind4, 0, imIndexToVector(y - disp, x + disp, h));
		//			valsG4(ind4)=dS(disp+1);
					cvmSet(valsG4, ind4, 0, cvmGet(dS, 0, disp));
		//		end;
				}
		//	end;
			}
		//end;
		}
		//%'done inds'
		//indsGx1=indsGx1(1:indx);
		//indsGx2=indsGx2(1:indx);
		//valsGx=valsGx(1:indx);
		//indsGy1=indsGy1(1:indy);
		//indsGy2=indsGy2(1:indy);
		//valsGy=valsGy(1:indy);
		//indsG31=indsG31(1:ind3);
		//indsG32=indsG32(1:ind3);
		//valsG3=valsG3(1:ind3);
		//indsG41=indsG41(1:ind4);
		//indsG42=indsG42(1:ind4);
		//valsG4=valsG4(1:ind4);

		//Gx=sparse(indsGx1,indsGx2,valsGx,imgSize,imgSize);
		sparse(indsGx1, indsGx2, valsGx, indx + 1, imgSize, imgSize, Gx);

		//Gy=sparse(indsGy1,indsGy2,valsGy,imgSize,imgSize);
		sparse(indsGy1, indsGy2, valsGy, indy + 1, imgSize, imgSize, Gy);

		//G3=sparse(indsG31,indsG32,valsG3,imgSize,imgSize);
		sparse(indsG31, indsG32, valsG3, ind3 + 1, imgSize, imgSize, G3);

		//G4=sparse(indsG41,indsG42,valsG4,imgSize,imgSize);
		sparse(indsG41, indsG42, valsG4, ind4 + 1, imgSize, imgSize, G4);

		cvReleaseMat(&dS);
		cvReleaseMat(&indsGx1);
		cvReleaseMat(&indsGx2);
		cvReleaseMat(&valsGx);
		cvReleaseMat(&indsGy1);
		cvReleaseMat(&indsGy2);
		cvReleaseMat(&valsGy);
		cvReleaseMat(&indsG31);
		cvReleaseMat(&indsG32);
		cvReleaseMat(&valsG3);
		cvReleaseMat(&indsG41);
		cvReleaseMat(&indsG42);
		cvReleaseMat(&valsG4);
	}

	//function [F,B]=solveFB(I,alpha)
	void solveFB(IplImage *src, IplImage *alpha, IplImage *&foreground, IplImage *&background) {
		//--PM(alpha);

		if (foreground)
			cvReleaseImage(&foreground);
		foreground = cvCloneImage(src);
		cvZero(foreground);

		if (background)
			cvReleaseImage(&background);
		background = cvCloneImage(src);
		cvZero(background);

		//[h,w,c]=size(I);
		int h = src->height;
		int w = src->width;
		int c = src->nChannels;
		
		//mask=(alpha>=0.02).*(alpha<=0.98);
		IplImage *mask = cvCreateImage(cvGetSize(src), 8, 1);	//	mask中1表示需要计算前景和背景区域的值，0表示不需要
		for (int ar = 0; ar < alpha->height; ar++)
			for (int ac = 0; ac < alpha->width; ac++)
				if ((0.02 * 255 <= CV_IMAGE_ELEM(alpha, BYTE, ar, ac)) && (CV_IMAGE_ELEM(alpha, BYTE, ar, ac) <= 0.98 * 255))
					CV_IMAGE_ELEM(mask, BYTE, ar, ac) = 1;
				else
					CV_IMAGE_ELEM(mask, BYTE, ar, ac) = 0;

		//--PM(mask);
		//SHOWTIME("fgbg, make mask");

		//[Gx,Gy,Gd1,Gd2]=getGMatByMask(w,h,mask);
		CvSparseMat *Gx = NULL;
		CvSparseMat *Gy = NULL;
		CvSparseMat *Gd1 = NULL;
		CvSparseMat *Gd2 = NULL;
		GetGMatByMask(w, h, mask, Gx, Gy, Gd1, Gd2);
		//SHOWTIME("fgbg, get Gmat");

		int sizes[2];

		//G=[Gx;Gy;Gd1;Gd2];
		CvSparseMat *G = NULL;
		sizes[0] = cvGetDimSize(Gx, 0) + cvGetDimSize(Gy, 0) + cvGetDimSize(Gd1, 0) + cvGetDimSize(Gd2, 0);
		sizes[1] = cvGetDimSize(Gx, 1);		
		G = cvCreateSparseMat(2, sizes, CV_64FC1);
		CvSparseMatIterator mat_iterator;	
		for (CvSparseNode* node = cvInitSparseMatIterator(Gx, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(Gx, node); 
			double val = *(double*)CV_NODE_VAL(Gx, node);		
			cvSetReal2D(G, idx[0], idx[1], val);
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(Gy, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(Gy, node); 
			double val = *(double*)CV_NODE_VAL(Gy, node);		
			cvSetReal2D(G, idx[0] + cvGetDimSize(Gx, 0), idx[1], val);
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(Gd1, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(Gd1, node); 
			double val = *(double*)CV_NODE_VAL(Gd1, node);		
			cvSetReal2D(G, idx[0] + cvGetDimSize(Gx, 0) + cvGetDimSize(Gy, 0), idx[1], val);
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(Gd2, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(Gd2, node); 
			double val = *(double*)CV_NODE_VAL(Gd2, node);		
			cvSetReal2D(G, idx[0] + cvGetDimSize(Gx, 0) + cvGetDimSize(Gy, 0) + cvGetDimSize(Gd1, 0), idx[1], val);
		}

		//SHOWTIME("Build G");
		//--PM(G, "G.txt");

		//Ga=G*alpha(:);
		CvMat *valpha = cvCreateMat(w * h, 1, CV_64FC1);
		MakeVector(alpha, valpha);

		CvMat *Ga = cvCreateMat(w * h * 4, 1, CV_64FC1);
		Mul(G, valpha, Ga);

		//--PM(Ga, "Ga.txt");
		//SHOWTIME("Build Ga");

		//wgf=abs(Ga).^0.5+0.003*repmat((1-alpha(:)),4,1);
		CvMat *wgf = cvCreateMat(w * h * 4, 1, CV_64FC1);
		cvZero(wgf);
		for (int i = 0; i < w * h * 4; i++)
			cvmSet(wgf, i, 0, sqrt(abs(cvmGet(Ga, i, 0))) + 0.003 * (1 - cvmGet(valpha, i % (w * h), 0)));

		//wgb=abs(Ga).^0.5+0.003*repmat(alpha(:),4,1);
		CvMat *wgb = cvCreateMat(w * h * 4, 1, CV_64FC1);
		cvZero(wgb);
		for (int i = 0; i < w * h * 4; i++)
			cvmSet(wgb, i, 0, sqrt(abs(cvmGet(Ga, i, 0))) + 0.003 * cvmGet(valpha, i % (w * h), 0));

		//wf=(alpha(:)>0.98)*100+0.03*alpha(:).*(alpha(:)>0.7)+0.01*(alpha(:)<0.02);
		CvMat *wf = cvCreateMat(w * h, 1, CV_64FC1);
		cvZero(wf);
		for (int i = 0; i < w * h; i++)
			cvmSet(wf, i, 0, 
				(cvmGet(valpha, i, 0) > 0.98) * 100 + 
				0.03 * cvmGet(valpha, i, 0) * (cvmGet(valpha, i, 0) > 0.7) +
				0.01 * (cvmGet(valpha, i, 0) < 0.02));

		//wb=(alpha(:)<0.02)*100+0.03*(1-alpha(:)).*(alpha(:)<0.3)+0.01*(alpha(:)>0.98);
		CvMat *wb = cvCreateMat(w * h, 1, CV_64FC1);
		cvZero(wb);
		for (int i = 0; i < w * h; i++)
			cvmSet(wb, i, 0, 
			(cvmGet(valpha, i, 0) < 0.02) * 100 + 
			0.03 * (1 - cvmGet(valpha, i, 0)) * (cvmGet(valpha, i, 0) < 0.3) +
			0.01 * (cvmGet(valpha, i, 0) > 0.98));

		int s_gbf[2] = {w * h * 4, w * h * 4};
		int s_bf[2] = {w * h, w * h};

		//wgf=spdiags(wgf(:),0,length(wgf),length(wgf));	
		CvSparseMat *s_wgf = cvCreateSparseMat(2, s_gbf, CV_64FC1);
		spdiags0(wgf, w * h * 4, s_wgf);

		//wgb=spdiags(wgb(:),0,length(wgb),length(wgb)); 
		CvSparseMat *s_wgb = cvCreateSparseMat(2, s_gbf, CV_64FC1);
		spdiags0(wgb, w * h * 4, s_wgb);

		//wf=spdiags(wf(:),0,length(wf),length(wf));
		CvSparseMat *s_wf = cvCreateSparseMat(2, s_bf, CV_64FC1);
		spdiags0(wf, w * h, s_wf);

		//wb=spdiags(wb(:),0,length(wb),length(wb)); 
		CvSparseMat *s_wb = cvCreateSparseMat(2, s_bf, CV_64FC1);
		spdiags0(wb, w * h, s_wb);

		//SHOWTIME("Build w");

		//--PM(s_wgf, "wgf.txt");
		//--PM(s_wgb, "wgb.txt");
		//--PM(s_wf, "wf.txt");
		//--PM(s_wb, "wb.txt");

		int dimensionsmall[2] = {4 * w * h, w * h};
		int dimensionmiddle1[2] = {2 * w * h, 2 * w * h};
		int dimensionmiddle2[2] = {1 * w * h, 2 * w * h};
		int dimensionmiddle3[2] = {8 * w * h, 2 * w * h};
		int dimensionlarge[2] = {11 * w * h, 2 * w * h};
		int dimensionlargeT[2] = {2 * w * h, 11 * w * h};
		int dimensionlargedenominator[2] = {2 * w * h, 2 * w * h};
		CvSparseMat *Ag00 = cvCreateSparseMat(2, dimensionsmall, CV_64FC1);
		CvSparseMat *Ag11 = cvCreateSparseMat(2, dimensionsmall, CV_64FC1);

		Mul(s_wgf, G, Ag00);
		Mul(s_wgb, G, Ag11);

		CvSparseMat *Ai = cvCreateSparseMat(2, dimensionmiddle1, CV_64FC1);
		CombineDiag(s_wf, s_wb, Ai);
		
		CvSparseMat *As = cvCreateSparseMat(2, dimensionmiddle2, CV_64FC1);
		for (int i = 0; i < w * h; i++) {
			cvSetReal2D(As, i, i, cvmGet(valpha, i, 0));
			cvSetReal2D(As, i, i + w * h, 1 - cvmGet(valpha, i, 0));
		}

		CvSparseMat *Ag = cvCreateSparseMat(2, dimensionmiddle3, CV_64FC1);
		CombineDiag(Ag00, Ag11, Ag);

		CvSparseMat *A = cvCreateSparseMat(2, dimensionlarge, CV_64FC1);	
		CombineRow(Ai, As, Ag, A);

		CvSparseMat *AT = cvCreateSparseMat(2, dimensionlargeT, CV_64FC1);	
		Transpose(A, AT);

		CvSparseMat *denominator = cvCreateSparseMat(2, dimensionlargedenominator, CV_64FC1);
		Mul(AT, A, denominator);

		//--PM(denominator, "deno.txt");
		//SHOWTIME("Build denominator");

		CvMat *bi = cvCreateMat(2 * w * h, 1, CV_64FC1);
		CvMat *bs = NULL;
		CvMat *bg = cvCreateMat(8 * w * h, 1, CV_64FC1);
		cvZero(bg);

		CvMat *tI = cvCreateMat(src->height, src->width, CV_64FC1);
		CvMat *vtI = cvCreateMat(w * h, 1, CV_64FC1);

		CvMat *bi0 = cvCreateMat(w * h, 1, CV_64FC1);
		CvMat *bi1 = cvCreateMat(w * h, 1, CV_64FC1);
		
		CvMat *b = cvCreateMat(11 * w * h, 1, CV_64FC1);

		CvMat *numerator = cvCreateMat(2 * w * h, 1, CV_64FC1);

		CvMat *x = cvCreateMat(11 * w * h, 1, CV_64FC1);

		//cout << "pre over" << endl;
		//for t=1:c
		for (int t = 0; t < c; t++ ) {
		//	tI=I(:,:,t);
			GetMatFromImage(src, t, tI);
			MakeVector(tI, vtI);

		//	Ag=[wgf*G,sparse(size(G,1),size(G,2));sparse(size(G,1),size(G,2)),wgb*G];
		//  bg=zeros(size(Ag,1),1);
		//  Ai=[wf,sparse(w*h,w*h);sparse(w*h,w*h),wb];
		//  bi=[wf*tI(:).*(alpha(:)>0.02);wb*tI(:).*(alpha(:)<0.98)];
			cvZero(bi);

			Mul(s_wf, vtI, bi0);
			for (int i = 0; i < w * h; i++)
				if (cvmGet(valpha, i, 0) > 0.02)
					cvmSet(bi, i, 0, cvmGet(bi0, i, 0));

			Mul(s_wb, vtI, bi1);
			for (int i = 0; i < w * h; i++)
				if (cvmGet(valpha, i, 0) < 0.98)
					cvmSet(bi, i + w * h, 0, cvmGet(bi1, i, 0));

		//  As=[spdiags(alpha(:),0,w*h,w*h),spdiags(1-alpha(:),0,w*h,w*h)];
		//  bs=tI(:);
			if (bs) cvReleaseMat(&bs);
			bs = cvCloneMat(vtI);

		//  A=[Ai;As;Ag];
		//  b=[bi;bs;bg];
			CombineRow(bi, bs, bg, b);

		//  x=(A'*A)\(A'*b);
			Mul(AT, b, numerator);

			//--if (t == 0)
			//--	PM(numerator, "numer0.txt");
			//--else if (t == 1)
			//--	PM(numerator, "numer1.txt");
			//--else 
			//--	PM(numerator, "numer2.txt");

			//SolveSymSparseMat(denominator, numerator, x, 2, 1e-6);
			KunXuSolver(denominator, numerator, x);

		//  F(:,:,t)=reshape(x(1:w*h),h,w);
			for (int i = 0; i < h; i++)
				for (int j = 0; j < w; j++) {
					CV_IMAGE_ELEM(foreground, BYTE, i, j * 3 + t) = Round(cvmGet(x, j * h + i, 0) * 255);
				}

		//  B(:,:,t)=reshape(x(w*h+1:end),h,w);
			for (int i = 0; i < h; i++)
				for (int j = 0; j < w; j++) {
					CV_IMAGE_ELEM(background, BYTE, i, j * 3 + t) = Round(cvmGet(x, j * h + i + w * h, 0) * 255);
				}
			//SHOWTIME(t);
		//end  
		}

		cvReleaseSparseMat(&Gx);
		cvReleaseSparseMat(&Gy);
		cvReleaseSparseMat(&Gd1);
		cvReleaseSparseMat(&Gd2);
		cvReleaseSparseMat(&G);
		cvReleaseMat(&valpha);
		cvReleaseMat(&Ga);
		cvReleaseMat(&wgf);
		cvReleaseMat(&wgb);
		cvReleaseMat(&wf);
		cvReleaseMat(&wb);
		cvReleaseSparseMat(&s_wgf);
		cvReleaseSparseMat(&s_wgb);
		cvReleaseSparseMat(&s_wf);
		cvReleaseSparseMat(&s_wb);
		cvReleaseSparseMat(&Ag00);
		cvReleaseSparseMat(&Ag11);
		cvReleaseSparseMat(&Ai);
		cvReleaseSparseMat(&As);
		cvReleaseSparseMat(&Ag);
		cvReleaseSparseMat(&A);
		cvReleaseSparseMat(&AT);
		cvReleaseSparseMat(&denominator);
		cvReleaseMat(&bi);
		cvReleaseMat(&bs);
		cvReleaseMat(&bg);
		cvReleaseMat(&tI);
		cvReleaseMat(&vtI);
		cvReleaseMat(&bi0);
		cvReleaseMat(&bi1);
		cvReleaseMat(&b);
		cvReleaseMat(&numerator);
		cvReleaseMat(&x);
		cvReleaseImage(&mask);
	}

	void GetGMatByMaskX(int w, int h, IplImage *mask, CvSparseMat *&Gx, CvSparseMat *&Gy, CvSparseMat *&G3, CvSparseMat *&G4) {
		int imgSize = w * h;
		int dims[2] = {imgSize, imgSize};
		if (Gx) cvReleaseSparseMat(&Gx);
		if (Gy) cvReleaseSparseMat(&Gy);
		if (G3) cvReleaseSparseMat(&G3);
		if (G4) cvReleaseSparseMat(&G4);
		Gx = cvCreateSparseMat(2, dims, CV_64FC1);
		Gy = cvCreateSparseMat(2, dims, CV_64FC1);
		G3 = cvCreateSparseMat(2, dims, CV_64FC1);
		G4 = cvCreateSparseMat(2, dims, CV_64FC1);

		for (int x = 0; x < w - 1; x++) {
		//	for y=1:h,
			for (int y = 0; y < h; y++) {
		//		if ((~mask(y,x))&(~mask(y,x+1)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y, x + 1))) {
		//			continue
					continue;
		//		end  
				}
		//		for disp=0:filtSizeS,
		//			indx=indx+1;
		//			indsGx1(indx)=imIndexToVect(y,x,h);				
		//			indsGx2(indx)=imIndexToVect(y,x+disp,h);
		//			valsGx(indx)=dS(disp+1);
		//		end
				int tr = imIndexToVector(y, x + 0, h);
				int tc = imIndexToVector(y, x + 1, h);
				cvSetReal2D(Gx, tr, tr, 1);
				cvSetReal2D(Gx, tr, tc, -1);
		//	end
			}
		//end
		}
		
		//for x=1:w,	
		for (int x = 0; x < w; x++) {		
		//	for y=1:h-1,
			for (int y = 0; y < h - 1; y++) {
		//		if ((~mask(y,x))&(~mask(y+1,x)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y + 1, x))) {
		//			continue
					continue;
		//		end
				}
		//		for disp=0:filtSizeS,
		//			indy=indy+1;
		//			indsGy1(indy)=imIndexToVect(y,x,h);
		//			indsGy2(indy)=imIndexToVect(y+disp,x,h);
		//			valsGy(indy)=dS(disp+1);
		//		end;
				int tr = imIndexToVector(y + 0, x, h);
				int tc = imIndexToVector(y + 1, x, h);
				cvSetReal2D(Gy, tr, tr, 1);
				cvSetReal2D(Gy, tr, tc, -1);
		//	end;
			}
		//end
		}

		//for x=1:w-1,	
		for (int x = 0; x < w - 1; x++) {		
		//	for y=1:h-1,
			for (int y = 0; y < h - 1; y++) {
		//		if ((~mask(y,x))&(~mask(y+1,x+1)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y + 1, x + 1))) {
		//			continue
					continue;
		//		end
				}
		//		for disp=0:filtSizeS,
		//			ind3=ind3+1;
		//			indsG31(ind3)=imIndexToVect(y,x,h);
		//			indsG32(indy)=imIndexToVect(y+disp,x+disp,h);
		//			valsG3(ind3)=dS(disp+1);
		//		end
				int tr = imIndexToVector(y + 0, x + 0, h);
				int tc = imIndexToVector(y + 1, x + 1, h);
				cvSetReal2D(G3, tr, tr, 1);
				cvSetReal2D(G3, tr, tc, -1);
		//	end
			}			
		//end
		}

		//for x=1:w-1,	
		for (int x = 0; x < w - 1; x++) {
		//	for y=2:h,
			for (int y = 1; y < h; y++) {
		//		if ((~mask(y,x))&(~mask(y-1,x+1)))
				if ((!CV_IMAGE_ELEM(mask, BYTE, y, x)) & (!CV_IMAGE_ELEM(mask, BYTE, y - 1, x + 1))) {
		//			continue
					continue;
		//		end
				}
		//		for disp=0:filtSizeS,   
		//			ind4=ind4+1;
		//			indsG41(ind4)=imIndexToVect(y,x,h);
		//			indsG42(ind4)=imIndexToVect(y-disp,x+disp,h);
		//			valsG4(ind4)=dS(disp+1);
		//		end;
				int tr = imIndexToVector(y - 0, x + 0, h);
				int tc = imIndexToVector(y - 1, x + 1, h);
				cvSetReal2D(G4, tr, tr, 1);
				cvSetReal2D(G4, tr, tc, -1);
		//	end;
			}
		//end;
		}
	}

	void getNonZeroX(CvSparseMat *Gx, CvSparseMat *Gy, CvSparseMat *G3, CvSparseMat *G4, IplImage *alpha, int *nonZeroIndex, int *&allNonZero, int &totalNonZero) {
		int maxNonZero = TotalNonZero(Gx) * 10;
		int *all = new int[maxNonZero];
		int pall = 0;

		CvSparseMatIterator mat_iterator;
		for (CvSparseNode* node = cvInitSparseMatIterator(Gx, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(Gx, node); 
			if (idx[0] != idx[1]) {
				all[pall++] = idx[0];
				all[pall++] = idx[1];
			}
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(Gy, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(Gy, node); 
			if (idx[0] != idx[1]) {
				all[pall++] = idx[0];
				all[pall++] = idx[1];
			}
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(G3, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(G3, node); 
			if (idx[0] != idx[1]) {
				all[pall++] = idx[0];
				all[pall++] = idx[1];
			}
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(G4, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(G4, node); 
			if (idx[0] != idx[1]) {
				all[pall++] = idx[0];
				all[pall++] = idx[1];
			}
		}

		for (int r = 0; r < alpha->height; r++) {
			for (int c = 0; c < alpha->width; c++) {
				if ((CV_IMAGE_ELEM(alpha, double, r, c) != 0) && (CV_IMAGE_ELEM(alpha, double, r, c) != 1)) {
					all[pall++] = r + c * alpha->height;					
				}
			}
		}
		qsort(all, pall, sizeof(int), qsortcompare);

		totalNonZero = 0;
		for (int i = 0; i < pall; i++) {
			if ((i == 0) || ((i != 0) && (all[i] != all[i - 1]))) {
				totalNonZero++;
			}
		}
		allNonZero = new int[totalNonZero];
		totalNonZero = 0;
		for (int i = 0; i < pall; i++) {
			if ((i == 0) || ((i != 0) && (all[i] != all[i - 1]))) {
				allNonZero[totalNonZero++] = all[i];
			}
		}

		for (int i = 0; i < totalNonZero; i++) {
			nonZeroIndex[allNonZero[i]] = i;
		}
	}

	void removeNonZeroX(CvSparseMat *&G, int *nonZeroIndex, int totalNonZero) {
		int dims[2] = {totalNonZero, totalNonZero};
		CvSparseMat *tg = cvCreateSparseMat(2, dims, CV_64FC1);
		CvSparseMatIterator mat_iterator;
		for (CvSparseNode* node = cvInitSparseMatIterator(G, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(G, node);
			double val = *(double*)CV_NODE_VAL(G, node);
			cvSetReal2D(tg, nonZeroIndex[idx[0]], nonZeroIndex[idx[1]], val);
		}
		cvReleaseSparseMat(&G);
		G = tg;
	}

	void buildWGX(CvSparseMat *G, CvMat *alpha, bool isfg, CvSparseMat *elem) {
		int total = cvGetSize(alpha).height;
		int dims[2] = {total, total};

		CvMat *ta = cvCloneMat(alpha);
		Mul(G, alpha, ta);
		SqrtAbs(ta);
		for (int i = 0; i < total; i++) {
			double v = 0;
			if (isfg) {
				v = cvmGet(ta, i, 0) + 0.003 * (1 - cvmGet(alpha, i, 0));
			}
			else {
				v = cvmGet(ta, i, 0) + 0.003 * cvmGet(alpha, i, 0);
			}
			cvmSet(ta, i, 0, v);
		}

		CvSparseMat *tg1 = cvCreateSparseMat(2, dims, CV_64FC1);
		Transpose(G, tg1);

		CvSparseMat *tg2 = cvCreateSparseMat(2, dims, CV_64FC1);
		for (int i = 0; i < total; i++) {
			cvSetReal2D(tg2, i, i, cvmGet(ta, i, 0) * cvmGet(ta, i, 0));
		}

		CvSparseMat *tg3 = cvCloneSparseMat(G);

		CvSparseMat *tmpelem = cvCreateSparseMat(2, dims, CV_64FC1);
		Mul(tg1, tg2, tmpelem);
		Mul(tmpelem, tg3, elem);
		cvReleaseSparseMat(&tmpelem);
		cvReleaseSparseMat(&tg1);
		cvReleaseSparseMat(&tg2);
		cvReleaseSparseMat(&tg3);
		cvReleaseMat(&ta);
	}

	void solveFBX(IplImage *src, IplImage *alpha, IplImage *&foreground, IplImage *&background) {
		if (foreground)
			cvReleaseImage(&foreground);
		foreground = cvCloneImage(src);
		cvZero(foreground);

		if (background)
			cvReleaseImage(&background);
		background = cvCloneImage(src);
		cvZero(background);

		int h = src->height;
		int w = src->width;
		int c = src->nChannels;
		int imgSize = h * w;

		IplImage *mask = cvCreateImage(cvGetSize(src), 8, 1);	//	mask中1表示需要计算前景和背景区域的值，0表示不需要
		for (int ar = 0; ar < alpha->height; ar++)
			for (int ac = 0; ac < alpha->width; ac++)
				if ((0.02 <= CV_IMAGE_ELEM(alpha, double, ar, ac)) && (CV_IMAGE_ELEM(alpha, double, ar, ac) <= 0.98))
					CV_IMAGE_ELEM(mask, BYTE, ar, ac) = 1;
				else
					CV_IMAGE_ELEM(mask, BYTE, ar, ac) = 0;
		
		CvSparseMat *Gx = NULL;	
		CvSparseMat *Gy = NULL;	
		CvSparseMat *G3 = NULL;
		CvSparseMat *G4 = NULL;
		GetGMatByMaskX(w, h, mask, Gx, Gy, G3, G4);	
		//SHOWTIME("Get G");
		
		int totalNonZero;
		int *allNonZero = NULL;
		int *nonZeroIndex = new int[h * w];
		getNonZeroX(Gx, Gy, G3, G4, alpha, nonZeroIndex, allNonZero, totalNonZero);

		int nonZeroDims[2] = {totalNonZero, totalNonZero};
		//SHOWTIME("Get getNonZeroX");

		removeNonZeroX(Gx, nonZeroIndex, totalNonZero);
		removeNonZeroX(Gy, nonZeroIndex, totalNonZero);
		removeNonZeroX(G3, nonZeroIndex, totalNonZero);
		removeNonZeroX(G4, nonZeroIndex, totalNonZero);

		CvMat *valpha = cvCreateMat(totalNonZero, 1, CV_64FC1);
		for (int i = 0; i < totalNonZero; i++) {
			int tr = allNonZero[i] % h;
			int tc = allNonZero[i] / h;			
			cvmSet(valpha, i, 0, CV_IMAGE_ELEM(alpha, double, tr, tc));
		}	
		double *wf = new double[totalNonZero];
		double *wb = new double[totalNonZero];

		for (int i = 0; i < totalNonZero; i++) {
			double ta = cvmGet(valpha, i, 0);
			double tf = (ta > 0.98) * 100 + 0.03 * ta * (ta > 0.7) + 0.01 * (ta < 0.02);
			double tb = (ta < 0.02) * 100 + 0.03 * (1 - ta) * (ta < 0.3) + 0.01 * (ta > 0.98);
			wf[i] = tf;
			wb[i] = tb;
		}

		CvSparseMat *Afg[4];
		for (int i = 0; i < 4; i++) { 
			Afg[i] = cvCreateSparseMat(2, nonZeroDims, CV_64FC1);	
		}		
		CvSparseMat *Abg[4];
		for (int i = 0; i < 4; i++) { 
			Abg[i] = cvCreateSparseMat(2, nonZeroDims, CV_64FC1);	
		}		

		buildWGX(Gx, valpha, true, Afg[0]);		
		buildWGX(Gy, valpha, true, Afg[1]);
		buildWGX(G3, valpha, true, Afg[2]);
		buildWGX(G4, valpha, true, Afg[3]);

		buildWGX(Gx, valpha, false, Abg[0]);
		buildWGX(Gy, valpha, false, Abg[1]);
		buildWGX(G3, valpha, false, Abg[2]);
		buildWGX(G4, valpha, false, Abg[3]); 

		int nonZeroDimsFB[2] = {totalNonZero * 2, totalNonZero * 2};
		CvSparseMat *A = cvCreateSparseMat(2, nonZeroDimsFB, CV_64FC1);

		for (int i = 0; i < 4; i++) {
			CvSparseMatIterator mat_iterator;
			for (CvSparseNode* node = cvInitSparseMatIterator(Afg[i], &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator))
			{
				int* idx = CV_NODE_IDX(Afg[i], node); 
				double val = *(double*)CV_NODE_VAL(Afg[i], node);		

				double v = cvGetReal2D(A, idx[0], idx[1]) + val;				
				cvSetReal2D(A, idx[0], idx[1], v);
			}
		}
		for (int i = 0; i < 4; i++) {
			CvSparseMatIterator mat_iterator;
			for (CvSparseNode* node = cvInitSparseMatIterator(Abg[i], &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator))
			{
				int* idx = CV_NODE_IDX(Abg[i], node); 
				double val = *(double*)CV_NODE_VAL(Abg[i], node);		

				double v = cvGetReal2D(A, idx[0] + totalNonZero, idx[1] + totalNonZero) + val;				
				cvSetReal2D(A, idx[0] + totalNonZero, idx[1] + totalNonZero, v);
			}
		}		

		//	完成A的左上块，主要是前景
		for (int i = 0; i < totalNonZero; i++) {
			double v = cvGetReal2D(A, i, i);
			v += wf[i] * wf[i] + cvmGet(valpha, i, 0) * cvmGet(valpha, i, 0);
			cvSetReal2D(A, i, i, v);
		}
		//	完成A的右下块，主要是背景
		for (int i = 0; i < totalNonZero; i++) {
			double v = cvGetReal2D(A, i + totalNonZero, i + totalNonZero);
			v += wb[i] * wb[i] + (1 - cvmGet(valpha, i, 0)) * (1 - cvmGet(valpha, i, 0));
			cvSetReal2D(A, i + totalNonZero, i + totalNonZero, v);
		}
		//	完成A的左下块和右上块
		for (int i = 0; i < totalNonZero; i++) {
			double v = cvmGet(valpha, i, 0) * (1 - cvmGet(valpha, i, 0));
			cvSetReal2D(A, i, i + totalNonZero, v);
			cvSetReal2D(A, i + totalNonZero, i, v);
		}

		CvMat *b = cvCreateMat(2 * totalNonZero, 1, CV_64FC1);
		CvMat *x = cvCreateMat(2 * totalNonZero, 1, CV_64FC1);

		for (int t = 0; t < 3; t++) {
			for (int i = 0; i < totalNonZero; i++) {
				int tr = allNonZero[i] % h;
				int tc = allNonZero[i] / h;
				double tI = CV_IMAGE_ELEM(src, BYTE, tr, tc * 3 + t) / 255.0;
				double va = cvmGet(valpha, i, 0);
				double v = 0;

				//	完成B的上块，主要是前景
				v = wf[i] * wf[i] * tI * (va > 0.02) + va * tI;
				cvmSet(b, i, 0, v);
				//	完成B的下块，主要是背景
				v = wb[i] * wb[i] * tI * (va < 0.98) + (1 - va) * tI;
				cvmSet(b, i + totalNonZero, 0, v);
			}		

			cvZero(x);
			KunXuSolver(A, b, x);

			for (int i = 0; i < totalNonZero; i++) {
				int tr = allNonZero[i] % h;
				int tc = allNonZero[i] / h;
				double v = 0;

				v = cvmGet(x, i, 0);
				CV_IMAGE_ELEM(foreground, BYTE, tr, tc * 3 + t) = Round(v * 255);
				v = cvmGet(x, i + totalNonZero, 0);
				CV_IMAGE_ELEM(background, BYTE, tr, tc * 3 + t) = Round(v * 255);
			}
		}

		int pall = 0;
		for (int c = 0; c < w; c++) {
			for (int r = 0; r < h; r++) {
				int idx = imIndexToVector(r, c, h);
				if (idx == allNonZero[pall]) {
					pall++;
				}		
				else {
					if (CV_IMAGE_ELEM(alpha, double, r, c) == 1.0) {
						CV_IMAGE_ELEM(foreground, BYTE, r, c * 3 + 0) = CV_IMAGE_ELEM(src, BYTE, r, c * 3 + 0);
						CV_IMAGE_ELEM(foreground, BYTE, r, c * 3 + 1) = CV_IMAGE_ELEM(src, BYTE, r, c * 3 + 1);
						CV_IMAGE_ELEM(foreground, BYTE, r, c * 3 + 2) = CV_IMAGE_ELEM(src, BYTE, r, c * 3 + 2);
					}
					else {
						CV_IMAGE_ELEM(background, BYTE, r, c * 3 + 0) = CV_IMAGE_ELEM(src, BYTE, r, c * 3 + 0);
						CV_IMAGE_ELEM(background, BYTE, r, c * 3 + 1) = CV_IMAGE_ELEM(src, BYTE, r, c * 3 + 1);
						CV_IMAGE_ELEM(background, BYTE, r, c * 3 + 2) = CV_IMAGE_ELEM(src, BYTE, r, c * 3 + 2);
					}
				}
			}
		}

		cvReleaseSparseMat(&A);
		cvReleaseMat(&b);
		cvReleaseMat(&x);
		for (int i = 0; i < 4; i++) {
			cvReleaseSparseMat(&(Afg[i]));
			cvReleaseSparseMat(&(Abg[i]));
		}

		cvReleaseImage(&mask);
		cvReleaseSparseMat(&Gx);
		cvReleaseSparseMat(&Gy);
		cvReleaseSparseMat(&G3);
		cvReleaseSparseMat(&G4);
		delete[] allNonZero;
		delete[] nonZeroIndex;
	}
}