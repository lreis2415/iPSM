#include "EnvDataset.h"

namespace solim {
	EnvDataset::EnvDataset()
		: LayerRef(nullptr), CellSize(-9999.), CellSizeY(-9999.), XMin(-9999.), XMax(-9999.),
		YMin(-9999.), YMax(-9999.), XStart(0), YStart(0), TotalX(0), TotalY(0), CalcArea(0) {
		LayerNames.clear();
		LayerNames.shrink_to_fit();
	}

    EnvDataset::EnvDataset(vector<string> &envLayerFilenames, vector<string> &datatypes){
		LayerNames.clear();
		LayerNames.shrink_to_fit();
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
			LayerNames.push_back(layername);
        }
        ReadinLayers(envLayerFilenames, datatypes, LayerNames, 1);
    }

	EnvDataset::EnvDataset(vector<string>& envLayerFilenames, vector<string>& datatypes, vector<string>& layernames, double ramEfficent)
		: LayerRef(nullptr), CellSize(-9999.), CellSizeY(-9999.), XMin(-9999.), XMax(-9999.),
		YMin(-9999.), YMax(-9999.), XStart(0), YStart(0), TotalX(0), TotalY(0), CalcArea(0) {
		LayerNames = layernames;
		LayerNames.shrink_to_fit();
		ReadinLayers(envLayerFilenames, datatypes, layernames, ramEfficent);
	}

	EnvDataset::~EnvDataset() {
	}


	void EnvDataset::RemoveAllLayers() {
		Layers.clear();
	}

	void EnvDataset::ReadinLayers(vector<string>& envLayerFilenames, const vector<string>& datatypes, vector<string>& layernames, double ramEfficent) {
		if (envLayerFilenames.empty() || datatypes.empty()) {
			// Print some error information and return.
			return;
		}
		int layerNum = int(envLayerFilenames.size());
		if (layerNum != int(datatypes.size())) {
			// Print some error information and return.
			return;
		}
		// Step 1. Read the header information of the first environment layer (as reference for comparison) using tiffIO
		EnvLayer *firstLayer = new EnvLayer(0, LayerNames[0], envLayerFilenames[0].c_str(), getDatatypeFromString(datatypes[0]));
		LayerRef = new BaseIO(firstLayer->baseRef);
		AddLayer(firstLayer);
		LayerNames.push_back(layernames[0]);

		TotalX = LayerRef->getXSize();
		TotalY = LayerRef->getYSize();
		CellSize = LayerRef->getDx();//LayerRef->getDxA();
		CellSizeY = LayerRef->getDy(); //LayerRef->getDyA();// Assuming dx==dy
		NoDataValue = LayerRef->getNoDataValue();

		// Read tiff data into partitions and blocks

		if (ramEfficent > 0.9999)
			LayerRef->blockNull();
		else
			LayerRef->blockInit(ramEfficent / double(layerNum));
		// Get the size of current block;
		XSize = LayerRef->getBlockX();
		YSize = LayerRef->getBlockY();
		//LayerRef->localToGlobal(0, 0, XStart, YStart);	// get the position of the current partition

		// get the global coordinates
		XMin = LayerRef->getXMin();
		YMax = LayerRef->getYMax();
		XMax = XMin + CellSize * TotalX;
		YMin = YMax - CellSizeY * TotalY;

		// Step 3. Create EnvLayer objects using linearpart data
		for (int i = 1; i < layerNum; ++i) {
            EnvLayer *newLayer = new EnvLayer(i, LayerNames[i], envLayerFilenames[i].c_str(), getDatatypeFromString(datatypes[i]), LayerRef);
			if (!LayerRef->compareIO(newLayer->baseRef)) {
                cout << "Warning: File needs to be reprojected: " << envLayerFilenames[i] << endl;
                vector<string> nameparts;
                ParseStr(envLayerFilenames[i],'.',nameparts);
                string resampleFile = "";
                for(size_t k = 0; k < nameparts.size() - 1; k++){
                    resampleFile +=nameparts[k];
                }
                if(nameparts.size()>1)
                    resampleFile = resampleFile+"_resampleForSoLIM."+nameparts[nameparts.size()-1];
                else
                    resampleFile = envLayerFilenames[i]+"_resampleForSoLIM.tif";
                bool success = Resample(envLayerFilenames[i],resampleFile, firstLayer);
                delete newLayer;
                if(success){
                    newLayer = new EnvLayer(i, LayerNames[i], resampleFile, getDatatypeFromString(datatypes[i]), LayerRef);
                    AddLayer(newLayer);
					LayerNames.push_back(layernames[i]);
                } else return;
			}
			else {
				AddLayer(newLayer);
				LayerNames.push_back(layernames[i]);
			}
		}
		LayerSize = Layers.size();
	}

	EnvUnit* EnvDataset::GetEnvUnit(const int row, const int col) {
		// receive global col and row number
		EnvUnit *e = new EnvUnit();
		e->Loc->Row = row;
		e->Loc->Col = col;
		e->Loc->X = col * CellSize + XMin;
		e->Loc->Y = YMax - row * CellSize;
		int numRows = 1;
		int numCols = 1;
		for (int i = 0; i < Layers.size(); ++i) {
			float *value = new float;
			Layers.at(i)->baseRef->read(e->Loc->Col, e->Loc->Row, numRows, numCols, value);
			e->AddEnvValue(Layers.at(i)->LayerName, *value, Layers.at(i)->DataType);
		}
		return e;
	}

	EnvUnit* EnvDataset::GetEnvUnit(const double x, const double y) {
		if (x<XMin || x>XMax || y<YMin || y>YMax) return nullptr;
		EnvUnit *e = new EnvUnit();
		e->Loc->X = x;
		e->Loc->Y = y;
		e->Loc->Row = int((YMax - y) / CellSize);
		e->Loc->Col = int((x - XMin) / CellSize);
		int numRows = 1;
		int numCols = 1;
		for (int i = 0; i < Layers.size(); ++i) {
			float *value = new float;
			*value = (float)this->NoDataValue;
			Layers.at(i)->baseRef->read(e->Loc->Col, e->Loc->Row, numRows, numCols, value);
			e->AddEnvValue(Layers.at(i)->LayerName, *value, Layers.at(i)->DataType);
		}
		return e;
	}

	void EnvDataset::GetEnvUnitValues(int valPos, float*envVals) {
		
		// receive global col and row number
		if (!envVals) envVals = new float[LayerSize];
		for (int i = 0; i < Layers.size(); ++i) {
			envVals[i] = Layers.at(i)->EnvData[valPos];
		}
	}

	EnvLayer *EnvDataset::getDEM() {
		for (auto it = Layers.begin(); it != Layers.end(); ++it) {
			string name = (*it)->LayerName;
			for (int i = 0; i < name.length(); ++i) {
				toupper(name[i]);
			}
			if (name == "DEM" || name == "ELEVATION") {
				return (*it);
			}
		}
		return nullptr;
	}
	void EnvDataset::Writeout(string filename, float* EnvData, int blockRank) {
		int localx = 0;
		int localy = 0;
		int globalx, globaly;
		LayerRef->localToGlobal(blockRank, localx, localy, globalx, globaly);
		int nx = LayerRef->getBlockX();
		int ny = LayerRef->getBlockY();
		if (blockRank == (LayerRef->getBlockSize() - 1)) {
			ny = LayerRef->getYSize() - blockRank * LayerRef->getBlockY();
		}
		if (blockRank == 0) LayerRef->writeInit();
		if (EnvData != nullptr)
			LayerRef->write(globalx, globaly, ny, nx, EnvData, filename);
	}


	bool EnvDataset::Resample(string& filename, string&resampled_fn, EnvLayer* refLayer) {
		//resample
		GDALAllRegister();
		GDALDataset* origin_ds = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
		if (origin_ds == NULL) {
			cout << "Error opening file " << filename << endl;
			return false;
		}
		GDALRasterBand* origin_band = origin_ds->GetRasterBand(1);
		GDALDataset* ref_ds = (GDALDataset*)GDALOpen(refLayer->FileName.c_str(), GA_ReadOnly);
		if (ref_ds == NULL) {
			cout << "Error opening file " << filename << endl;
			GDALClose(origin_ds);
			return false;
		}
		GDALRasterBand* ref_band = ref_ds->GetRasterBand(1);
		char* cFileName = new char[resampled_fn.length() + 1];
		strcpy(cFileName, resampled_fn.c_str());
		GDALDataType eBDataType = origin_band->GetRasterDataType();
		GDALDataset* dest_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(cFileName,
			TotalX, TotalY, 1,
			eBDataType, NULL);
		GDALRasterBand* dst_band = dest_ds->GetRasterBand(1);
		dst_band->SetColorInterpretation(origin_band->GetColorInterpretation());
		int success = -1;
		double noDataValue = NODATA;
		noDataValue = origin_band->GetNoDataValue(&success);
		if (!success)
		{
			GDALClose(origin_ds);
			GDALClose(dest_ds);
			GDALClose(ref_ds);
			return false;
		}
		dst_band->SetNoDataValue(noDataValue);
		GDALColorTable* colorTable = origin_band->GetColorTable();
		if (colorTable)
		{
			dst_band->SetColorTable(colorTable);
		}
		dest_ds->SetProjection(ref_ds->GetProjectionRef());
		double ref_adfGeoTransform[6];
		ref_ds->GetGeoTransform(ref_adfGeoTransform);
		dest_ds->SetGeoTransform(ref_adfGeoTransform);

		CPLErr result = GDALReprojectImage(origin_ds, origin_ds->GetProjectionRef(), dest_ds, ref_ds->GetProjectionRef(), GRA_NearestNeighbour, 0.0, 0.0, NULL, NULL, NULL);
		if (result == CE_Failure) {
			cout << "Error resampling file " << filename << endl;
			GDALClose(ref_ds);
			GDALClose(origin_ds);
			GDALClose(dest_ds);
			return false;
		}
		GDALClose(ref_ds);
		GDALClose(origin_ds);
		GDALClose(dest_ds);
		return true;

	}
}
