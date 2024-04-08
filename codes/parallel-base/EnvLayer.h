/*!
 * @brief Environmental layer
 * @todo
 * @version 1.0
 * @revision  17-11-21 zhanglei - initial version
 *            17-11-27 lj       - code revision and format
 */
#ifndef ENVLAYER_HPP_
#define ENVLAYER_HPP_

#include <string>
#include <vector>
//#include <gdal_priv.h>
#include "rasterLayer.h"
#include "DataTypeEnum.h"
#define VERY_SMALL 0.0001
#define NODATA -9999
// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::string;
using std::vector;
using GPRO::RasterLayer;
using GPRO::CoordBR;

namespace solim {
class EnvLayer {
public:
    int LayerId;
	RasterLayer<float> *baseRef;
    string LayerName;
	string FileName;
	float* EnvData;
	float* upperBorder;
	float* lowerBorder;
	int borderRowNum;
    // Adding MembershipData to use save membershipValue
    float* MembershipData;
    DataTypeEnum DataType;
    double Data_Max;
    double Data_Min;
    double Data_Range;
	double Data_Mean;
	double Data_StdDev; // standard deviation
    float NoDataValue;
    double GeoTransform[6];
	double CellSize;
    int XSize;
    int YSize;

public:
	EnvLayer();
	EnvLayer(const int layerId, string layerName, const string& filename, const DataTypeEnum dataType);
    ~EnvLayer();
    
    void read();
	void read(CoordBR *subWorkBR);
	bool compareIO(const EnvLayer *layer);
	void WriteoutMembership(const string& filename);
    void Writeout(string filename);
	void CalcStat();
	void CopyFrame(EnvLayer *lyr);
};
}

#endif
