#include "stdafx.h"
#include "Global.h"

namespace alphamatting {

	int compelemrowfirst(const sparselem *p1, const sparselem *p2) {
		if (p1->r < p2->r)
			return -1;
		else if (p1->r == p2->r && p1->c < p2->c)
			return -1;
		else if (p1->r == p2->r && p1->c == p2->c)
			return 0;
		else
			return 1;
	}

	int compelemcolfirst(const sparselem *p1, const sparselem *p2) {
		if (p1->c < p2->c)
			return -1;
		else if (p1->c == p2->c && p1->r < p2->r)
			return -1;
		else if (p1->c == p2->c && p1->r == p2->r)
			return 0;
		else
			return 1;
	}

	void PM(CvMat *m, char *filename) {
		if (filename) {
			FILE *f = fopen(filename, "w");
			CvSize size = cvGetSize(m);
			char buf[50];
			for (int i = 0; i < size.height; i++) {
				fprintf(f, "{");
				for (int j = 0; j < size.width; j++) {
					sprintf(buf, "%f", cvmGet(m, i, j));
					fprintf(f, "%15s, ", buf);
				}
				fprintf(f, "},");
				fprintf(f, "\n");
			}
			fclose(f);
		}
		else {
			CvSize size = cvGetSize(m);
			for (int r = 0; r < size.height; r++) {
				for (int c = 0; c < size.width; c++)
					cout << cvmGet(m, r, c) << ' ';
				cout << endl;
			}
		}
	}

	void PM(CvSparseMat *m, char *filename) {
		int total = GetSparseNumberOfElements(m);
		sparselem *elems = new sparselem[total];
		GetSparseElements(m, elems);
		qsort(elems, total, sizeof(sparselem), (int (*)(const void *, const void *))compelemcolfirst);

		if (filename) {
			FILE *f = fopen(filename, "w");
			//char buf[200];
			for (int i = 0; i < total; i++)
				fprintf(f, "%0.4d %0.4d %f\n", elems[i].r, elems[i].c, elems[i].v);
			fclose(f);
		}
		else {
			for (int i = 0; i < total; i++)
				cout << elems[i].r << ' ' << elems[i].c << ' ' << elems[i].v << endl;
		}

		delete[] elems;
	}

	void PM(IplImage *m, char *filename) {
		if (filename) {
			FILE *f = fopen(filename, "w");
			CvSize size = cvGetSize(m);
			char buf[50];
			for (int i = 0; i < size.height; i++) {
				fprintf(f, "{");
				for (int j = 0; j < size.width; j++) {
					sprintf(buf, "%f", CV_IMAGE_ELEM(m, BYTE, i, j));
					fprintf(f, "%15s, ", buf);
				}
				fprintf(f, "},");
				fprintf(f, "\n");
			}
			fclose(f);
		}
		else {
			CvSize size = cvGetSize(m);
			for (int r = 0; r < size.height; r++) {
				for (int c = 0; c < size.width; c++)
					cout << (float)CV_IMAGE_ELEM(m, BYTE, r, c) << ' ';
				cout << endl;
			}
		}
	}

}