#pragma once
class ParamLoader
{
public:
	ParamLoader(void);
	~ParamLoader(void);

	static cvi* LoadParamFromDir(string paramDir, string imgPrefix);
	static void SaveParamToDir(string paramDir, string imgPrefix, cvi* param);
	static void ShowParamInDir(string paramDir, string imgPrefix, cvi* param);

	static void SaveCviCh(int width, int height, cvi* param, int ch, string fileDir);
	static cvi* VislzCviCh(int width, int height, cvi* param, int ch, double minV, double maxV);

};

