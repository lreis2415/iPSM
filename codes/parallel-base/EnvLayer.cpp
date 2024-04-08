#include "EnvLayer.h"

//#include "ogrsf_frmts.h"  // for ogr
//#include "gdal_alg.h"     // for GDALPolygonize
//#include "cpl_conv.h"     // for CPLMalloc()

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::string;

namespace solim {
EnvLayer::EnvLayer()
    : LayerId(-1), DataType(CONTINUOUS), baseRef(nullptr), 
      Data_Max(NODATA), Data_Min(NODATA),Data_Range(NODATA), Data_Mean(NODATA),
	  Data_StdDev(NODATA), NoDataValue(NODATA),EnvData(nullptr),
      XSize(-1), YSize(-1) {
    for (int i = 0; i < 6; i++) GeoTransform[i] = NODATA;
}

EnvLayer::EnvLayer(const int layerId, string layerName, const string& filename, const DataTypeEnum dataType):
    LayerId(layerId), LayerName(layerName), DataType(dataType), FileName(filename){
	baseRef = new RasterLayer<float>(layerName);
	baseRef->readGlobalFileSerial(filename.c_str());
	double fileinfo[15];
	if(0==GetRank()){
		XSize = baseRef->_pMetaData->column; fileinfo[0] = XSize;
		YSize = baseRef->_pMetaData->row; fileinfo[1] = YSize;
		CellSize = baseRef->_pMetaData->cellSize; fileinfo[2] = CellSize;
		NoDataValue = baseRef->_pMetaData->noData; fileinfo[3] = NoDataValue;
		Data_Max = baseRef->dataMax; fileinfo[4] = Data_Max;
		Data_Mean = baseRef->dataMean; fileinfo[5] = Data_Mean;
		Data_Min = baseRef->dataMin; fileinfo[6] = Data_Min;
		Data_Range = baseRef->dataRange; fileinfo[7] = Data_Range;
		Data_StdDev = baseRef->dataStdDev; fileinfo[8] = Data_StdDev;
		for (int i = 0; i < 6; i++) { GeoTransform[i] = baseRef->_pMetaData->pTransform[i]; fileinfo[9+i] = GeoTransform[i];}
		MPI_Bcast(fileinfo, 15, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		MPI_Bcast(fileinfo, 15, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		XSize = fileinfo[0];
		YSize = fileinfo[1];
		CellSize = fileinfo[2];
		NoDataValue = fileinfo[3];
		Data_Max = fileinfo[4];
		Data_Mean = fileinfo[5];
		Data_Min = fileinfo[6];
		Data_Range = fileinfo[7];
		Data_StdDev = fileinfo[8];
		for (int i = 0; i < 6; i++) { GeoTransform[i] = fileinfo[9+i];} 
	}
	EnvData = nullptr;
}


EnvLayer::~EnvLayer(void) {
	delete baseRef;
}

void EnvLayer::read() {
	baseRef->readFile(FileName.c_str(),GPRO::ROWWISE_DCMP);
	EnvData = baseRef->cellSpace()->_matrix;
}

void EnvLayer::read(CoordBR *subWorkBR) {
	baseRef->readFile(FileName.c_str(),*subWorkBR, GPRO::ROWWISE_DCMP);
	EnvData = baseRef->cellSpace()->_matrix;
}

bool EnvLayer::compareIO(const EnvLayer *layer) {
	// compare the coordinates extent of the two layers
	if (XSize != layer->XSize || YSize != layer->YSize) {
		cout << "Columns or Rows do not match" << endl;
		return false;
	}
	if (abs(CellSize - layer->CellSize) > VERY_SMALL) {
		cout << "cellsize do not match" << endl;
		return false;
	}
	double xLeftEdge = GeoTransform[0];
	double yTopEdge = GeoTransform[3];
	if (abs(xLeftEdge - layer->GeoTransform[0]) > CellSize || abs(yTopEdge - layer->GeoTransform[3]) > CellSize) {
		cout << "Extent coordinate does not match exactly" << endl;
		return false;
	}
	return true;
}

void EnvLayer::CalcStat() {
	//todo
}

void EnvLayer::WriteoutMembership(const string& filename) {
    //todo
}

void EnvLayer::Writeout(string filename) {
	baseRef->writeFile(filename.c_str());
}

void EnvLayer::CopyFrame(EnvLayer *lyr) {
	CellSize = lyr->CellSize;
	XSize = lyr->XSize;
	YSize = lyr->YSize;
	baseRef = new RasterLayer<float>;
    baseRef->copyLayerInfo(*(lyr->baseRef));
	EnvData = baseRef->cellSpace()->_matrix;
	//upperBorder = new float[XSize];
	//lowerBorder = new float[XSize];
	borderRowNum = lyr->borderRowNum;
	DataType = CONTINUOUS;
	NoDataValue = lyr->NoDataValue;
	for (int i = 0; i<6; i++) GeoTransform[i] = lyr->GeoTransform[i];
}

};
