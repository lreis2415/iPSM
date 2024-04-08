#include "EnvDataset.h"

//#include "TauDEM_IO/createpart.h"

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::transform;

namespace solim {
EnvDataset::EnvDataset()
    : LayerRef(nullptr), CellSize(-9999.), CellSizeY(-9999.), XMin(-9999.), XMax(-9999.),
      YMin(-9999.), YMax(-9999.), XStart(0), YStart(0), TotalX(0), TotalY(0), CalcArea(0), LayerSize(0) {
}

EnvDataset::EnvDataset(vector<string>& envLayerFilenames, vector<string>& datatypes)
    : LayerRef(nullptr), CellSize(-9999.), CellSizeY(-9999.), XMin(-9999.), XMax(-9999.),
      YMin(-9999.), YMax(-9999.), XStart(0), YStart(0), TotalX(0), TotalY(0), CalcArea(0), LayerSize(0) {
	vector<string> layernames;
	for (size_t i = 0; i < envLayerFilenames.size(); i++) {
		string layername = "";
		string filename = envLayerFilenames[i];
		if (!filename.empty()) {
			std::size_t first = filename.find_last_of('/');
			if (first == std::string::npos) {
				first = filename.find_last_of('\\');
			}
			std::size_t end = filename.find_last_of('.');
			if (end == std::string::npos) end = filename.size();
			layername = filename.substr(first + 1, end - first - 1).c_str();
		}
		layernames.push_back(layername);
	}
    ReadinLayers(envLayerFilenames, datatypes, layernames);
}

EnvDataset::EnvDataset(vector<string> &envLayerFilenames, vector<string> &datatypes, vector<string>& layernames)
	: LayerRef(nullptr), CellSize(-9999.), CellSizeY(-9999.), XMin(-9999.), XMax(-9999.),
	YMin(-9999.), YMax(-9999.), XStart(0), YStart(0), TotalX(0), TotalY(0), CalcArea(0), LayerSize(0) {
	ReadinLayers(envLayerFilenames, datatypes, layernames);
}

EnvDataset::~EnvDataset() {
}


void EnvDataset::RemoveAllLayers() {
    Layers.clear();
}

void EnvDataset::AddLayer(EnvLayer *layer){
	if (LayerSize == 0) {
		Layers.push_back(layer);
		LayerNames.push_back(layer->LayerName);
	}
	else {
		if (!Layers.at(0)->compareIO(layer)) return;
		Layers.push_back(layer);
		LayerNames.push_back(layer->LayerName);
	}
	LayerSize++;
}
    
void EnvDataset::ReadinLayers(vector<string>& envLayerFilenames, const vector<string>& datatypes, const vector<string>& layernames) {
    /* Use the parallel framework of TauDEM to read raster data */
    if (envLayerFilenames.empty() || datatypes.empty()) {
		cout << "No files read";
        return;
    }
    int layerNum = int(envLayerFilenames.size());
    if (layerNum != int(datatypes.size())) {
        // Print some error information and return.
        return;
    }
    // Step 1. Read the header information of the environment layer using RasterLayer
	EnvLayer *firstLayer = new EnvLayer(0, layernames[0], envLayerFilenames[0], getDatatypeFromString(datatypes[0]));
	LayerRef = firstLayer->baseRef;
	TotalX = firstLayer->XSize;
	TotalY = firstLayer->YSize;
	CellSize = firstLayer->CellSize;
	CellSizeY = firstLayer->CellSize; // Assuming dx==dy
	NoDataValue = firstLayer->NoDataValue;
	//get the global coordinates
	XMin = firstLayer->GeoTransform[0];
	YMax = firstLayer->GeoTransform[3];
	XMax = XMin + CellSize * TotalX;
	YMin = YMax - CellSizeY * TotalY;

	// Step 2. Divide the raster into partitions based on PaRGO v2
	LoadRef = new ComputeLayer<float>("LayerRef");
	LoadRef->init(LayerRef, nullptr, 10);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	Transformation<float> transOper(0, 1, LoadRef);
	transOper.run();
	CoordBR subWorkBR;
	LoadRef->getCompuLoad(GPRO::ROWWISE_DCMP, process_nums, subWorkBR); // Decompose the spatial computational domain.
	AddLayer(firstLayer);
	for (int i = 1; i < layerNum; i++) {
		AddLayer(new EnvLayer(i, layernames[i], envLayerFilenames[i], getDatatypeFromString(datatypes[i])));
	}
	// Step 3. Read file
	for (int i = 0; i < LayerSize; i++) {
		Layers[i]->read(&subWorkBR);
	}
    // Get the size of current partition
	XSize = LayerRef->_pMetaData->_localdims.nCols();
    YSize = LayerRef->_pMetaData->_localdims.nRows();
	XStart = LayerRef->_pMetaData->_MBR.minICol();
	YStart = LayerRef->_pMetaData->_MBR.minIRow();
	LocalXMin = LocalXMin+ LayerRef->_pMetaData->_MBR.minICol() * CellSize;
	LocalYMax = LocalYMax - LayerRef->_pMetaData->_MBR.minIRow() * CellSizeY;
	LocalXMax = XMin + CellSize * XSize;
	LocalYMin = YMax - CellSizeY * YSize;
    
}


EnvUnit* EnvDataset::GetEnvUnit(const int row, const int col) {
	// receive global col and row number
	EnvUnit *e = new EnvUnit();
	e->Loc->Row = row;
	e->Loc->Col = col;
	e->Loc->X = col * CellSize + XMin;
	e->Loc->Y = YMax - row * CellSize;
	if (row < 0 || row > TotalY || col < 0 || col > TotalX) {
		return nullptr;
	}
#ifdef RasterLayer_H
	float *envVals = new float[LayerSize];
	float *tmpVals = new float[LayerSize];
	if (row >= YStart && row < YStart + YSize && col >= XStart &&col < XStart + XSize) {
	// get EnvUnit info and send messages to other ranks
		int valPos = XSize*(row - YStart) + (col - XStart);
		for (int i = 0; i < Layers.size(); ++i) {
			tmpVals[i] = Layers.at(i)->EnvData[valPos];
		}
	}
	else {
		for (int i = 0; i < Layers.size(); ++i) tmpVals[i] = 0;
	}
	MPI_Allreduce(tmpVals, envVals, LayerSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	for (int i = 0; i < Layers.size(); ++i) {
		e->AddEnvValue(LayerNames[i], envVals[i], Layers.at(i)->DataType);
	}
#endif // RasterLayer_H
	return e;
}

EnvUnit* EnvDataset::GetEnvUnit(const double x, const double y) {
    if (x < XMin || x > XMax || y < YMin || y > YMax) {
        return nullptr;
    }
    int irow = int((YMax - y) / CellSize);
    int icol = int((x - XMin) / CellSize);
    return GetEnvUnit(irow, icol);
}

void EnvDataset::GetEnvUnitValues(const int valPos, float*envVals) {
	if(!envVals) envVals = new float[LayerSize];
	for (int i = 0; i < LayerSize; i++) {
		envVals[i] = Layers.at(i)->EnvData[valPos];
	}
}
}
