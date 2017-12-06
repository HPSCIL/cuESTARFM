#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "gdalwarper.h"
#include <stdio.h>
#include"math.h"
#include"trans.h"
#define num_thread 256
#define num_block 128
int main()
{
	/*const char* tifFile1="/media/pxt/18AEAECFAEAEA52A/data/newdata/2009-249-flaash.dat";
    const char* tifFile2="/media/pxt/18AEAECFAEAEA52A/data/newdata/2009-329-flaash.dat";
	const char* modFile1="/media/pxt/18AEAECFAEAEA52A/data/newdata/MODO9A1.A2009249.dat";
	const char* modFile2="/media/pxt/18AEAECFAEAEA52A/data/newdata/MODO9A1.A2009329-0.dat";

	const char* modFile0="/media/pxt/18AEAECFAEAEA52A/data/newdata/MODO9A1.A2009289.dat";
	const char* out="/media/pxt/18AEAECFAEAEA52A/data/newdata/2009-289new-flaash.tif";*/
	//const char* modFile1="D:\\data\\newdata\\MODO9A1.A2009249.dat";
	//const char* tifFile1="D:\\data\\newdata\\2009-249-flaash.dat";
	//const char* modFile2="D:\\data\\newdata\\MODO9A1.A2009329-0.dat";
	//const char* tifFile2="D:\\data\\newdata\\2009-329-flaash.dat";
	//const char* modFile0="D:\\data\\newdata\\MODO9A1.A2009289.dat";
	//const char* out="D:\\data\\newdata\\ESTARFM.tif";
	const char* modFile1="D:\\cuda\\shikong\\软件\\测试数据\\M_2002_01_04.tif";
    const char* modFile2="D:\\cuda\\shikong\\软件\\测试数据\\M_2002_02_21.tif";
	const char* tifFile1="D:\\cuda\\shikong\\软件\\测试数据\\L_2002_01_04.tif";
	const char* tifFile2="D:\\cuda\\shikong\\软件\\测试数据\\L_2002_02_21.tif";

	const char* modFile0="D:\\cuda\\shikong\\软件\\测试数据\\M_2002_02_12.tif";
	const char* out="D:\\cuda\\shikong\\软件\\测试数据\\Llick_2001_11_14.tif";
	float uncertain=0.0028;
	int w=51;
	int class_num=6;
	Re_fusion2(tifFile1, modFile1,tifFile2,modFile2,modFile0,out,w,class_num,uncertain);
}