#include "ipsmOperator.h"

#include <iomanip>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::setw;

bool solim::ipsmOperator::PredictMap_Property() {
	int xSize = EDS->TotalX;
	int ySize = EDS->TotalY;
	int count = xSize * ySize;
	if (nullptr == Map_Prediction) Map_Prediction = new EnvLayer();
	if (nullptr == Map_Uncertainty) Map_Uncertainty = new EnvLayer();
	Map_Prediction->CopyFrame(EDS->Layers[0]);
	Map_Uncertainty->CopyFrame(EDS->Layers[0]);

	float* data_prediction = Map_Prediction->EnvData;
	float* data_uncertainty = Map_Uncertainty->EnvData;

	int envUnitCount = EDS->XSize*EDS->YSize;//allEnvUnits.size();
	int sampleCount = SampleEnvUnits.size();
	int segementCount = 1000;
	EnvUnit *se = nullptr;
#ifdef BLOCKIO
	for (int blockNum = 0; blockNum < EDS->LayerRef->getBlockSize(); blockNum++) {
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
				cout << "Completed " << setw(5) << (int(((i*1.0+ blockNum*envUnitCount) / count + 0.5 / segementCount)*segementCount)) / (segementCount / 100.0) << "%";
#else
				cout << "Completed " << setw(5) << (int((i*1.0 / count + 0.5 / segementCount)*segementCount)) / (segementCount / 100.0) << "%";
#endif
			}
			float *envVals = new float[EDS->LayerSize];
			EDS->GetEnvUnitValues(i, envVals);
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
				double envSimi = CalcSimi(se, envVals);
				if (envSimi > simi_max) {
					simi_max = envSimi;
				}
				if (envSimi > thred_envsimi) {
					sum1 += envSimi * se->SoilVariable;
					sum2 += envSimi;
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
		Map_Prediction->Writeout(map_predict_filename,blockNum);
		Map_Uncertainty->Writeout(map_uncer_filename,blockNum);
	}
#else
	Map_Prediction->Writeout(map_predict_filename);
	Map_Uncertainty->Writeout(map_uncer_filename);
#endif // BLOCKIO

	se = nullptr;


	return true;
}
