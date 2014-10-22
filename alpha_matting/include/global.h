#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>


using namespace std;

extern int tm;
#ifndef SHOWTIME
#define SHOWTIME(x)	cout << x << ": " << GetTickCount() - tm << endl; \
					tm = GetTickCount()
#endif

namespace alphamatting {
	struct sparselem {
		int r, c;
		double v;
	};

	void GetRect(IplImage *imggray, int row, int col, int radius, RECT *lprect);
	void GetRect(CvMat *src, RECT *lprect, CvMat *dst);
	void GetRect(IplImage *img, RECT *lprect, int width, int height, CvMat *dst);
	int Round(double v);
	int GetSparseNumberOfElements(CvSparseMat *M);
	void GetSparseElements(CvSparseMat *M, sparselem *elems);
	void spdiags0(CvMat *m, int size, CvSparseMat *s);
	void MakeVector(CvMat *m, CvMat *v);
	void MakeVector(IplImage *img, CvMat *v);
	void Mul(CvSparseMat *A, CvSparseMat *B, CvSparseMat *C);
	void Mul(CvSparseMat *A, CvMat *B, CvMat *C);
	void SqrtAbs(CvMat *m);
	void CombineDiag(CvSparseMat *m00, CvSparseMat *m11, CvSparseMat *m);
	void CombineRow(CvMat *m1, CvMat *m2, CvMat *m3, CvMat *m);
	void CombineRow(CvSparseMat *m1, CvSparseMat *m2, CvSparseMat *m3, CvSparseMat *m);void Transpose(CvSparseMat *m, CvSparseMat *mT);
	void GetMatFromImage(IplImage *img, int ch, CvMat *m);
	int TotalNonZero(CvSparseMat *m);	
	void KunXuSolver(CvSparseMat *A, CvMat *b, CvMat *x);

	int compelemrowfirst(const sparselem *p1, const sparselem *p2);
	int compelemcolfirst(const sparselem *p1, const sparselem *p2);
	void PM(CvMat *m, char *filename = NULL);
	void PM(CvSparseMat *m, char *filename = NULL);
	void PM(IplImage *m, char *filename = NULL);
	void PM01(IplImage *m, char *filename = NULL);
	void PMline(CvMat *m, char *filename = NULL);

	void MakeFBUmask(IplImage *maskrgb, IplImage *&res, int unknownsize = 5);
	void GetLaplacian(IplImage *img, IplImage *consts, int win_size, CvSparseMat *&A, int *&nonZeroIndex, int *&allNonZero);

	/***********************************************************
	*	solveAlpha			���ݱ�����alpha
	*	���룺
	*			src			ԭʼRGBͼ��
	*			trimap		�Ҷ�ͼ��ǰ������Ϊ255����������Ϊ0
	*	�����
	*			alpha		�����ĻҶ�ͼ��0~255ӳ�䵽0~1
	*			alpha01		0~1֮��ʵ���洢�ĻҶ�ͼ
	***********************************************************/
	void solveAlpha(IplImage *src, IplImage *trimap, IplImage *&alpha, IplImage *&alpha01);

	/***********************************************************
	*	solveFB				����ԭͼ��alpha���ǰ���ͱ�������ȫ����matlab�㷨��д����
	*	���룺
	*			src			ԭʼRGBͼ��
	*			alpha		�Ҷ�ͼ��ʾ��alphaֵ��0~255ӳ�䵽0~1
	*	�����
	*			foreground	������ǰ����RGBͼ
	*			background	�����ı�����RGBͼ
	***********************************************************/
	void solveFB(IplImage *src, IplImage *alpha, IplImage *&foreground, IplImage *&background);

	/***********************************************************
	*	solveFBX			����ԭͼ��alpha���ǰ���ͱ������Ż���matlab�㷨��д����
	*	���룺
	*			src			ԭʼRGBͼ��
	*			alpha		�Ҷ�ͼ��ʾ��alphaֵ��0~255ӳ�䵽0~1
	*	�����
	*			foreground	������ǰ����RGBͼ
	*			background	�����ı�����RGBͼ
	***********************************************************/
	void solveFBX(IplImage *src, IplImage *alpha, IplImage *&foreground, IplImage *&background);
}