#include "StdAfx.h"
#include "param_loader.h"


ParamLoader::ParamLoader(void)
{
}


ParamLoader::~ParamLoader(void)
{
}

cvi* ParamLoader::LoadParamFromDir(string paramDir, string imgPrefix)
{
	string fileDir[3];
	fileDir[0] = paramDir + imgPrefix + "_0_gain.txt";
	fileDir[1] = paramDir + imgPrefix + "_1_bais.txt";
	fileDir[2] = paramDir + imgPrefix + "_2_bais.txt";

	int w, h;
	ifstream fin1(fileDir[0].c_str());
	if(!fin1)
	{
		cout<<"LoadGTFile: "<<fileDir[0]<<" not exist! error!"<<endl;
		return 0;
	}
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

void ParamLoader::SaveParamToDir(string paramDir, string imgPrefix, cvi* param)
{
	int w = param->width, h = param->height;
	SaveCviCh(w, h, param, 0, paramDir + imgPrefix + "_0_gain.txt");
	SaveCviCh(w, h, param, 1, paramDir + imgPrefix + "_1_bais.txt");
	SaveCviCh(w, h, param, 2, paramDir + imgPrefix + "_2_bais.txt");
}


void ParamLoader::SaveCviCh(int width, int height, cvi* param, int ch, string fileDir)
{
	ofstream fout(fileDir.c_str());
	fout<<width<<" "<<height<<endl;
	doF(i, height)
	{
		//cout<<"\r"<<i<<" ";
		doF(j, width)
		{ 
			double v = cvg2(param, i, j).val[ch];
			if(fequal(_f v, 0.f)) v = 0;
			fout<<v<<" ";
		}
	}
	fout.close();
}

void ParamLoader::ShowParamInDir(string paramDir, string imgPrefix, cvi* param)
{
	doF(k, 3)
	{
		cvi* temp = VislzCviCh(param->width, param->height, param, k, ((k==0)?0:-100), ((k==0)?1:100));
		cvsi(paramDir + imgPrefix + "_" + toStr(k) + ((k==0)?"_gain":"_bais") + ".png", temp); cvri(temp);
	}
}

cvi* ParamLoader::VislzCviCh(int width, int height, cvi* param, int ch, double minV, double maxV)
{
	cvi* res = cvci81(width, height);
	doFcvi(res, i, j)
	{
		double v = cvg2(param, i, j).val[ch];
		v = (v - minV) / (maxV - minV) * 255;
		cvs20(res, i, j, v);
	}
	return res;
}