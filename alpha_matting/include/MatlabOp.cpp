#include "stdafx.h"
#include "Global.h"
#include "../../matrix/matrix.h"

//#pragma comment(lib, "libguide40.lib")

namespace alphamatting {

	void GetRect(IplImage *imggray, int row, int col, int radius, RECT *lprect) {
		lprect->top = ((row - radius) >= 0) ? (row - radius) : 0;
		lprect->bottom = ((row + radius) < imggray->height) ? (row + radius) : imggray->height - 1;
		lprect->left = ((col - radius) >= 0) ? (col - radius) : 0;
		lprect->right = ((col + radius) < imggray->width) ? (col + radius) : imggray->width - 1;
	}

	void GetRect(CvMat *src, RECT *lprect, CvMat *dst) {
		int pid = 0;	
		for (int c = lprect->left; c <= lprect->right; c++)
			for (int r = lprect->top; r <= lprect->bottom; r++) {			
				cvmSet(dst, pid, 0, cvmGet(src, r, c));
				pid++;
			}
	}

	void GetRect(IplImage *img, RECT *lprect, int width, int height, CvMat *dst) {
		int channels = img->nChannels;
		int pos = 0;
		for (int ch = 0; ch < channels; ch++) {
			for (int c = lprect->left; c <= lprect->right; c++)
				for (int r = lprect->top; r <= lprect->bottom; r++) {				
					int dr = pos % height;
					int dc = pos / height;
					double v = CV_IMAGE_ELEM(img, BYTE, r, c * 3 + 2 - ch) / 255.0;
					cvmSet(dst, dr, dc, v);
					pos++;
				}
		}	
	}

	int GetSparseNumberOfElements(CvSparseMat *M) {
		CvSparseMatIterator mat_iterator;

		int pid = 0;
		for (CvSparseNode* node = cvInitSparseMatIterator(M, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator))
			pid++;
		return pid;
	}

	void GetSparseElements(CvSparseMat *M, sparselem *elems) {
		CvSparseMatIterator mat_iterator;

		int pid = 0;
		for (CvSparseNode* node = cvInitSparseMatIterator(M, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator))
		{
			int* idx = CV_NODE_IDX(M, node); 
			double val = *(double*)CV_NODE_VAL(M, node);		
			elems[pid].r = idx[0];
			elems[pid].c = idx[1];
			elems[pid].v = val;
			pid++;
		}
	}

	void MakeVector(CvMat *m, CvMat *v) {
		int w = cvGetSize(m).width;
		int h = cvGetSize(m).height;

		for (int c = 0; c < w; c++)
			for (int r = 0; r < h; r++)
				cvmSet(v, c * h + r, 0, cvmGet(m, r, c));
	}

	void MakeVector(IplImage *img, CvMat *v) {
		for (int r = 0; r < img->height; r++)
			for (int c = 0; c < img->width; c++)
				cvmSet(v, c * img->height + r, 0, CV_IMAGE_ELEM(img, BYTE, r, c) / 255.0);
	}

	void Mul(CvSparseMat *A, CvSparseMat *B, CvSparseMat *C) {
		int dimension[2] = {cvGetDimSize(A, 0), cvGetDimSize(B, 1)};
		int mid = cvGetDimSize(A, 1);

		int Atotal = GetSparseNumberOfElements(A);
		int Btotal = GetSparseNumberOfElements(B);

		sparselem *Aelems = new sparselem[Atotal];
		sparselem *Belems = new sparselem[Btotal];

		GetSparseElements(A, Aelems);
		GetSparseElements(B, Belems);

		qsort(Aelems, Atotal, sizeof(sparselem), (int (*)(const void *, const void *))compelemcolfirst);
		qsort(Belems, Btotal, sizeof(sparselem), (int (*)(const void *, const void *))compelemrowfirst);

		int *Aentrys = new int[mid + 1];
		int *Bentrys = new int[mid + 1];

		int point = 0;
		for (int i = 0; i < mid + 1; i++)	Aentrys[i] = Atotal;
		for (int i = 0; i < Atotal; i++) {
			if (i == 0 || Aelems[i].c != Aelems[i - 1].c) {			
				Aentrys[Aelems[i].c] = i;
			}
		}
		for (int i = mid - 1; i > 0; i--)
			if (Aentrys[i] == Atotal)
				Aentrys[i] = Aentrys[i + 1];
			
		point = 0;
		for (int i = 0; i < mid + 1; i++)	Bentrys[i] = Btotal;
		for (int i = 0; i < Btotal; i++) {
			if (i == 0 || Belems[i].r != Belems[i - 1].r) {			
				Bentrys[Belems[i].r] = i;
			}
		}
		for (int i = mid - 1; i > 0; i--)
			if (Bentrys[i] == Btotal)
				Bentrys[i] = Bentrys[i + 1];

		for (int i = 0; i < mid; i++) {
			for (int j = Aentrys[i]; j < Aentrys[i + 1]; j++)
				for (int k = Bentrys[i]; k < Bentrys[i + 1]; k++) {
					double v = cvGetReal2D(C, Aelems[j].r, Belems[k].c);
					v += Aelems[j].v * Belems[k].v;
					cvSetReal2D(C, Aelems[j].r, Belems[k].c, v);
				}				
		}

		delete[] Aelems;
		delete[] Belems;
		delete[] Aentrys;
		delete[] Bentrys;
	}

	void SqrtAbs(CvMat *m) {
		CvSize s = cvGetSize(m);
		for (int r = 0; r < s.height; r++)
			for (int c = 0; c < s.width; c++) {
				cvmSet(m, r, c, sqrt(fabs(cvmGet(m, r, c))));
			}
	}

	void Mul(CvSparseMat *A, CvMat *B, CvMat *C) {
		cvZero(C);

		CvSparseMatIterator mat_iterator;

		for (CvSparseNode* node = cvInitSparseMatIterator(A, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator))
		{
			int* idx = CV_NODE_IDX(A, node); 
			double val = *(double*)CV_NODE_VAL(A, node);		

			double v = cvmGet(C, idx[0], 0);
			v += val * cvmGet(B, idx[1], 0);
			cvmSet(C, idx[0], 0, v);
		}
	}

	void spdiags0(CvMat *m, int size, CvSparseMat *s) {
		for (int i = 0; i < size; i++)
			cvSetReal2D(s, i, i, cvmGet(m, i, 0));
	}

	void GetMatFromImage(IplImage *img, int ch, CvMat *m) {
		for (int r = 0; r < img->height; r++)
			for (int c = 0; c < img->width; c++) {
				cvmSet(m, r, c, CV_IMAGE_ELEM(img, BYTE, r, c * img->nChannels + ch) / 255.0);
			}
	}

	void CombineDiag(CvSparseMat *m00, CvSparseMat *m11, CvSparseMat *m) {
		int h_offset = cvGetDimSize(m00, 0);
		int w_offset = cvGetDimSize(m00, 1);

		CvSparseMatIterator mat_iterator;
		for (CvSparseNode* node = cvInitSparseMatIterator(m00, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(m00, node); 
			double val = *(double*)CV_NODE_VAL(m00, node);		
			cvSetReal2D(m, idx[0], idx[1], val);
		}	
		for (CvSparseNode* node = cvInitSparseMatIterator(m11, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(m11, node); 
			double val = *(double*)CV_NODE_VAL(m11, node);		
			cvSetReal2D(m, idx[0] + h_offset, idx[1] + w_offset, val);
		}	
	}

	void CombineRow(CvMat *m1, CvMat *m2, CvMat *m3, CvMat *m) {
		int h1 = cvGetSize(m1).height;
		int h2 = cvGetSize(m2).height;
		int h3 = cvGetSize(m3).height;

		int w = cvGetSize(m1).width;

		for (int i = 0; i < h1; i++)
			for (int j = 0; j < w; j++)
				cvmSet(m, i, j, cvmGet(m1, i, j));
		for (int i = 0; i < h2; i++)
			for (int j = 0; j < w; j++)
				cvmSet(m, i + h1, j, cvmGet(m2, i, j));
		for (int i = 0; i < h3; i++)
			for (int j = 0; j < w; j++)
				cvmSet(m, i + h1 + h2, j, cvmGet(m3, i, j));
	}

	void CombineRow(CvSparseMat *m1, CvSparseMat *m2, CvSparseMat *m3, CvSparseMat *m) {
		int h1 = cvGetDimSize(m1, 0);
		int h2 = cvGetDimSize(m2, 0);

		CvSparseMatIterator mat_iterator;
		for (CvSparseNode* node = cvInitSparseMatIterator(m1, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(m1, node); 
			double val = *(double*)CV_NODE_VAL(m1, node);		
			cvSetReal2D(m, idx[0], idx[1], val);
		}	
		for (CvSparseNode* node = cvInitSparseMatIterator(m2, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(m2, node); 
			double val = *(double*)CV_NODE_VAL(m2, node);		
			cvSetReal2D(m, idx[0] + h1, idx[1], val);
		}
		for (CvSparseNode* node = cvInitSparseMatIterator(m3, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(m3, node); 
			double val = *(double*)CV_NODE_VAL(m3, node);		
			cvSetReal2D(m, idx[0] + h1 + h2, idx[1], val);
		}
	}

	void Transpose(CvSparseMat *m, CvSparseMat *mT) {
		CvSparseMatIterator mat_iterator;
		for (CvSparseNode* node = cvInitSparseMatIterator(m, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(m, node); 
			double val = *(double*)CV_NODE_VAL(m, node);				
			cvSetReal2D(mT, idx[1], idx[0], val);
		}
	}

	int TotalNonZero(CvSparseMat *m) {
		int total = 0;
		CvSparseMatIterator mat_iterator;
		for (CvSparseNode* node = cvInitSparseMatIterator(m, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			total++;
		}
		return total;
	}

	int Round(double v) {
		if (v > 255)
			return 255;
		else if (v < 0)
			return 0;
		else
			return (int)floor(v + 0.5);
	}

	void KunXuSolver(CvSparseMat *A, CvMat *b, CvMat *x) {
		int m_size = cvGetDimSize(b, 0);

		sparse::symmatrix m;
		m.create(m_size);

		CvSparseMatIterator mat_iterator;	
		for (CvSparseNode* node = cvInitSparseMatIterator(A, &mat_iterator); node != 0; node = cvGetNextSparseNode(&mat_iterator)) {
			int* idx = CV_NODE_IDX(A, node); 
			double val = *(double*)CV_NODE_VAL(A, node);
			if (idx[0] <= idx[1])
				m.add(idx[0], idx[1], val);
		}
		
		double *tb = new double[m_size];
		double *tx = new double[m_size];
		for (int i = 0; i < m_size; i++)
			tb[i] = cvmGet(b, i, 0);
		
		m.solve(tx, tb);
		for (int i = 0; i < m_size; i++)
			cvmSet(x, i, 0, tx[i]);
	}
}