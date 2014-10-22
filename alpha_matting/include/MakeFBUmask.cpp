#include "stdafx.h"
#include "Global.h"

#define C_UNKNOW	1
#define C_FG		2
#define C_BG		0
#define C_BLACK		0
#define C_WHITE		2

namespace alphamatting {

	int dirs[][2] = {{0, -1}, {1, 0}, {0, 1}, {-1, 0}, {-1, -1}, {-1, 1}, {1, 1}, {1, -1}};

	const BYTE colors[][3] = {
		{0, 0, 0}, 
		{128, 128, 128},
		{255, 255, 255},
		{0, 0, 255}, 
		{255, 0, 0}, 
		{0, 255, 0}};

	bool IsColor(IplImage *img, int x, int y, int cid) {
		if (CV_IMAGE_ELEM(img, BYTE, y, x * 3 + 0) == colors[cid][2] &&
			CV_IMAGE_ELEM(img, BYTE, y, x * 3 + 1) == colors[cid][1] &&
			CV_IMAGE_ELEM(img, BYTE, y, x * 3 + 2) == colors[cid][0])
			return true;
		else
			return false;	
	}

	void SetColor(IplImage *img, int x, int y, int cid) {
		CV_IMAGE_ELEM(img, BYTE, y, x * 3 + 0) = colors[cid][2];
		CV_IMAGE_ELEM(img, BYTE, y, x * 3 + 1) = colors[cid][1];
		CV_IMAGE_ELEM(img, BYTE, y, x * 3 + 2) = colors[cid][0];
	}


	void circleit(IplImage *img, int x, int y, int radius, int cid) {
		for (int r = - radius; r <= radius; r++) 
			if (0 <= y + r && y + r < img->height) 
				for (int c = - radius; c <= radius; c++) 
					if (0 <= x + c && x + c < img->width)
						if (r * r + c * c <= radius * radius)
						{
							SetColor(img, x + c, y + r, cid);
							//CV_IMAGE_ELEM(back, BYTE, y + r, (x + c) * 3 + 0) = 0;
							//CV_IMAGE_ELEM(back, BYTE, y + r, (x + c) * 3 + 1) = 255;
							//CV_IMAGE_ELEM(back, BYTE, y + r, (x + c) * 3 + 2) = 0;
						}
	}

	bool IsBorder(IplImage *img, int r, int c, int cid, int ds = 4) {
		if (IsColor(img, c, r, cid)) {
			for (int d = 0; d < ds; d++)
				if (0 <= r + dirs[d][0] && r + dirs[d][0] < img->height &&
					0 <= c + dirs[d][1] && c + dirs[d][1] < img->width)
					if (IsColor(img, c + dirs[d][1], r + dirs[d][0], cid) == false)
						return true;
		}
		return false;
	}

	int total = 0;

	void AddOneMask(IplImage *bdmask, int r, int c, vector<CvPoint> &border) {
		total++;

		if (CV_IMAGE_ELEM(bdmask, BYTE, r, c) == 0)
			return;

		border.push_back(cvPoint(c, r));
		CV_IMAGE_ELEM(bdmask, BYTE, r, c) = 0;

		while (true) {
			int d;
			for (d = 0; d < 8; d++)
				if (0 <= r + dirs[d][0] && r + dirs[d][0] < bdmask->height &&
					0 <= c + dirs[d][1] && c + dirs[d][1] < bdmask->width)
					if (CV_IMAGE_ELEM(bdmask, BYTE, r + dirs[d][0], c + dirs[d][1]) == 255) {
						break;
					}
					if (d == 8)
						break;
					r += dirs[d][0];
					c += dirs[d][1];
					border.push_back(cvPoint(c, r));
					CV_IMAGE_ELEM(bdmask, BYTE, r, c) = 0;
		}
		border.push_back(cvPoint(-1, -1));
	}

	void GetBorder(IplImage *img, int cid, vector<CvPoint> &border) {
		IplImage *bdmask = cvCreateImage(cvGetSize(img), 8, 1);
		cvZero(bdmask);

		int fr;// fc;
		for (fr = 0; fr < img->height; fr++) {
			for (int fc = 0; fc < img->width; fc++) {
				if (IsBorder(img, fr, fc, cid))
					CV_IMAGE_ELEM(bdmask, BYTE, fr, fc) = 255;
			}
		}

		for (int i = 0; i < bdmask->width; i++) 
			if (CV_IMAGE_ELEM(bdmask, BYTE, 0, i) == 255) 
				AddOneMask(bdmask, 0, i, border);
		for (int i = 0; i < bdmask->width; i++) 
			if (CV_IMAGE_ELEM(bdmask, BYTE, bdmask->height - 1, i) == 255) 
				AddOneMask(bdmask, bdmask->height - 1, i, border);
		for (int i = 0; i < bdmask->height; i++) 
			if (CV_IMAGE_ELEM(bdmask, BYTE, i, 0) == 255) 
				AddOneMask(bdmask, i, 0, border);
		for (int i = 0; i < bdmask->height; i++) 
			if (CV_IMAGE_ELEM(bdmask, BYTE, i, bdmask->width - 1) == 255) 
				AddOneMask(bdmask, i, bdmask->width - 1, border);

		for (fr = 0; fr < img->height; fr++) {
			for (int fc = 0; fc < img->width; fc++) {
				if (CV_IMAGE_ELEM(bdmask, BYTE, fr, fc) == 255)
					AddOneMask(bdmask, fr, fc, border);				
			}
		}	
		cvReleaseImage(&bdmask);
	}


	void GetBorder(IplImage *mask, vector<CvPoint> &borders) {
		GetBorder(mask, C_FG, borders);
	}

	void MakeFBUmask(IplImage *maskrgb, IplImage *&res, int unknownsize) {
		vector<CvPoint> borders;
		GetBorder(maskrgb, borders);
		if (res)
			cvReleaseImage(&res);
		res = cvCloneImage(maskrgb);
		
		for (int i = 0; i < (int)borders.size(); i++)
			if (borders[i].x != -1 && borders[i].y != -1)
				circleit(res, borders[i].x, borders[i].y, unknownsize, C_UNKNOW);
	}

}