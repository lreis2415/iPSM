#include "ipsmIDWOperator.h"

#include <iomanip>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::setw;

bool solim::ipsmIDWOperator::PredictMap_Property() {
	int xSize = EDS->TotalX;
	int ySize = EDS->TotalY;
	double Cellsize = EDS->CellSize;
	double halfcell = Cellsize / 2;
	double Xmin = EDS->XMin;
	double Ymin = EDS->YMin;
	int iRow, iCol;
	double iX, iY, sX, sY;
	int count = xSize * ySize;

	if (nullptr == Map_Prediction) Map_Prediction = new EnvLayer();
	if (nullptr == Map_Uncertainty) Map_Uncertainty = new EnvLayer();
	Map_Prediction->CopyFrame(EDS->Layers.at(0));
	Map_Uncertainty->CopyFrame(EDS->Layers.at(0));

	float* data_prediction = Map_Prediction->EnvData;
	float* data_uncertainty = Map_Uncertainty->EnvData;
	
	int envUnitCount = EDS->XSize*EDS->YSize;//allEnvUnits.size();
	int sampleCount = SampleEnvUnits.size();
	double* dis = new double[sampleCount];
	double* simi_dis = new double[sampleCount];
	int segementCount = 1000;
	EnvUnit *se = nullptr;

#ifdef BLOCKIO
	for (int blockNum = 0; blockNum < EDS->LayerRef->getBlockSize(); blockNum++) {
		cout << "block num:" << blockNum << " of " << EDS->LayerRef->getBlockSize() << endl;
		for (int i = 0; i < EDS->Layers.size(); ++i) {
			EDS->Layers.at(i)->ReadByBlock(blockNum);
		}
		envUnitCount = xSize*EDS->Layers.at(0)->GetYSizeByBlock(blockNum);
#endif // BLOCKIO
		for (int i = 0; i < envUnitCount; i++)
		{
			
			if (i % (count / segementCount) == 0)
			{
				cout << '\r';
#ifdef BLOCKIO
				cout << "Completed " << setw(5) << (int(((i*1.0 + blockNum*envUnitCount) / count + 0.5 / segementCount)*segementCount)) / (segementCount / 100.0) << "%";
#else
				cout << "Completed " << setw(5) << (int((i*1.0 / count + 0.5 / segementCount)*segementCount)) / (segementCount / 100.0) << "%";
#endif
			}
			float *envVals = new float[EDS->LayerSize];
			EDS->GetEnvUnitValues(i, envVals);
			iRow = i / xSize;
			iCol = i % xSize;
			iX = Xmin+iCol*Cellsize + halfcell;
			//iY = Ymin+(ySize-iRow-1)*Cellsize + halfcell;
			//cout << "loc iY:" << iRow << endl;
#ifdef BLOCKIO
			iY = Ymin + (ySize - blockNum*blockYsize-iRow - 1)*Cellsize + halfcell;
#else
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			iY = Ymin + (ySize - EDS->LayerRef->_pMetaData->_MBR.minIRow() - iRow - 1)*Cellsize + halfcell;		
			//cout << "glb iY:" << EDS->LayerRef->_pMetaData->_MBR.minIRow() << endl;
#endif
			bool isCal = true;
			for (int k = 0; k < EDS->LayerSize; k++) {
				if (envVals[k] - EDS->Layers.at(k)->NoDataValue < VERY_SMALL)
					isCal = false;
			}
			if(!isCal) {
				data_prediction[i] = EDS->NoDataValue;//.push_back(EDS->NoDataValue);
				data_uncertainty[i] = EDS->NoDataValue;//.push_back(EDS->NoDataValue);
				continue;
			}

			double simi_max = 0.0;
			double sum1 = 0;	// sum of soil property * environmental similarity
			double sum2 = 0;	// sum of environmental similarities

			int tmpCount = 0;
			for (int j = 0; j < sampleCount; j++)
			{
				se = SampleEnvUnits[j];
				sX=se->Loc->X;
				sY=se->Loc->Y;
				double envSimi = CalcSimi(se, envVals, EuclideanDistance);;
				if (envSimi > simi_max) {
					simi_max = envSimi;
				}
				//cout << "envSimi" << envSimi << endl;
				if (envSimi > thred_envsimi) {
					
					double x_distance = sX - iX;
					double y_distance = sY - iY;
					dis[j] = sqrt(x_distance*x_distance + y_distance*y_distance);
					if ((fabs(x_distance) < halfcell) && (fabs(y_distance) < halfcell))
						simi_dis[j] = 1;		
					else
						simi_dis[j] = pow(dis[j], -r);
					sum1 += envSimi * se->SoilVariable*simi_dis[j];
					sum2 += envSimi*simi_dis[j] ;
					tmpCount++;
				}
			}

			if (tmpCount > 0) {
				data_prediction[i] = 1.0 * sum1 / sum2;//.push_back( 1.0 * sum1 / sum2 );
				data_uncertainty[i] = 1 - simi_max;//.push_back( 1 - simi_max );
			}
			else {
				data_prediction[i] = val_cannot_pred; //.push_back( val_cannot_pred );
				data_uncertainty[i] = 1 - simi_max;//.push_back( 1 - simi_max );
			}
		}
		cout << '\r';
		cout << "Completed " << setw(5) << 100.0 << "%";

#ifdef BLOCKIO
		Map_Prediction->Writeout(map_predict_filename, blockNum);
		Map_Uncertainty->Writeout(map_uncer_filename, blockNum);
	}
#else
	Map_Prediction->Writeout(map_predict_filename);
	Map_Uncertainty->Writeout(map_uncer_filename);
#endif // BLOCKIO

	se = nullptr;


	return true;
}
