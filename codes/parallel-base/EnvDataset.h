/*!
 * @brief Environmental dataset.
 * @todo 
 * @version 1.0
 * @revision  17-11-21 zhanglei - initial version
 *             17-11-27 lj       - code revision and format
 */
#ifndef ENVDATASET_HPP_
#define ENVDATASET_HPP_
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "computeLayer.h"
#include "transformation.h"
#include "EnvLayer.h"
#include "EnvUnit.h"
// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::transform;
using GPRO::ComputeLayer;
using GPRO::Transformation;

namespace solim {
class EnvDataset {
public:
    RasterLayer<float>* LayerRef;
	ComputeLayer<float>* LoadRef;
    vector<EnvLayer *> Layers;
    vector<string> LayerNames;
    // Basic Settings
    double CellSize;
    double CellSizeY;
    double XMin;
    double XMax;
    double YMin;
    double YMax;
	double LocalXMin, LocalXMax, LocalYMin, LocalYMax;
    double NoDataValue;
	int LayerSize;
    int XSize; // cols of current partition
    int YSize; // rows of current partition
    int XStart; // global coordinate of current partition's first cell's position
    int YStart;
    long TotalX; // global size
    long TotalY;
    int CalcArea;
	int rank, process_nums; // MPI params

public:
    EnvDataset();

    EnvDataset(vector<string> &envLayerFilenames, vector<string> &datatypes);
	EnvDataset(vector<string> &envLayerFilenames, vector<string> &datatypes, vector<string>& layernames);

    ~EnvDataset();

	void AddLayer(EnvLayer* layer);

    void RemoveAllLayers();

    void ReadinLayers(vector<string> &envLayerFilenames, const vector<string> &datatypes, const vector<string>& layernames);

	EnvLayer* getDEM();

    EnvUnit* GetEnvUnit(const int row, const int col);

    EnvUnit* GetEnvUnit(const double x, const double y);

	void GetEnvUnitValues(const int valPos, float *envVals);

	void Writeout(string filename, float*EnvData);

private:
    EnvDataset(const EnvDataset&);

    EnvDataset& operator=(const EnvDataset&);
};
}

#endif
