#include "ipsmAIDWOperator.h"

#include <iomanip>
#include <math.h>

#define PI 2*asin(1)
//using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::setw;

double AdaptiveR(vector<double> *a,  double expR, double *maxdis, double *mindis/*, int nbr*/){
	int count = (*a).size();
	vector<double> t = *a;
	double medianDis, r;
	double c;
	float con = 10;
	int nbr =5;
	if (nbr == 10)
		nbr = count / 2;
	for (int i = 0; i < count - 1; i++){
		for (int j = i+1; j < count; j++){
			if (t[i]>t[j]){
				c = t[i];
				t[i] = t[j];
				t[j] = c;
			}
		}
	}
	*maxdis = t[count * 9 / 10] / expR;
	*mindis = t[count / 10] / expR;
	//double maxDis = t[count * 9 / 10];
	//double minDis = t[count / 10];
	double dis=0;
	if (!(count % 2))
		medianDis = (t[count / 2] + t[count / 2 - 1]) / 2;
	else
		medianDis = t[(count - 1) / 2];
	if (count < nbr)
		nbr = count;
	for (int k = 0; k < nbr; k++){
		dis += t[k];
	}
	double obsR = dis / nbr;
	if (nbr == 5)
		obsR = medianDis;
	double R = obsR / expR;

	//for (int k = 0; k < count; k++)
	//	cout << "sample dis:" << t[k] << endl;
	
	
	return R;
}
double AdaptiveR2(const vector<double> a, int num, float con=10.5) {
	int count = a.size();
	vector<double> t = a;
	double medianDis, r;
	double c;
	//con = 1.;
	for (int i = 0; i < count - 1; i++) {
		for (int j = i + 1; j < count; j++) {
			if (t[i]>t[j]) {
				c = t[i];
				t[i] = t[j];
				t[j] = c;
			}
		}
	}
	//for (int k = 0; k < count; k++)
	//	cout << "sample dis:" << t[k] << endl;
	if (!(num % 2))
		medianDis = (t[num / 2] + t[num / 2 - 1]) / 2;
	else
		medianDis = t[(num - 1) / 2];
	//cout << "mediandis:" << medianDis << endl;
	r = log(2) / log(medianDis) * con;
	//r = log(2) / log(medianDis);
	//cout << "adaR:" << r << endl;
	return r;
}
double Trimember(double r, double maxR, double minR, double l1, double l5){
	double R;
	if (r >= maxR)
		R = 1;
	if (r < minR)
		R = 0;
	R = 0.5 + 0.5*sin((PI / (maxR - minR)*(r - minR)) - PI / 2);
	double  l2, l3, l4;
	double m;
	//l1 = 0.5;
	l2 = l1 + (l5 - l1)*0.25;
	l3 = l1 + (l5 - l1)*0.5;
	l4 = l1 + (l5 - l1)*0.75;
	//l5 = 0.75;
	if (R <= 0.)
		m = 1 * l1;
	else if (R <= 0.25)
		m = l2*(R - 0.) * 5 + l1*(1 - (R - 0.) * 5);
	else if (R <= 0.5)
		m = l3*(R - 0.25) * 5 + l2*(1 - (R - 0.25) * 5);
	else if (R <= 0.75)
		m = l4*(R - 0.5) * 5 + l3*(1 - (R - 0.5) * 5);
	else if (R <= 1.0)
		m = l5*(R - 0.75) * 5 + l4*(1 - (R - 0.75) * 5);
	else
		m = 1 * l5;
	return m;
}


bool solim::ipsmAIDWOperator::PredictMap_Property() {
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
	Map_Prediction->CopyFrame(EDS->Layers[0]);
	Map_Uncertainty->CopyFrame(EDS->Layers[0]);

	float* data_prediction = Map_Prediction->EnvData;
	float* data_uncertainty = Map_Uncertainty->EnvData;


	double simi_dis;
	//double maxR = 74, minR = 8;
	int envUnitCount = EDS->XSize*EDS->YSize;//allEnvUnits.size();
	int sampleCount = SampleEnvUnits.size();
	double* dis = new double[sampleCount];
	int segementCount = 1000;
	EnvUnit *se = nullptr;
	double expR = 1 / (2 * pow(sampleCount / (envUnitCount*Cellsize*Cellsize), 0.5));

	int c = 0;
#ifdef BLOCKIO
	for (int blockNum = 0; blockNum < EDS->LayerRef->getBlockSize(); blockNum++) {
		//cout << "block num:" << blockNum << " of " << EDS->LayerRef->getBlockSize() << endl;
		for (int i = 0; i < EDS->Layers.size(); ++i) {
			EDS->Layers.at(i)->ReadByBlock(blockNum);
		}
		envUnitCount = xSize*EDS->Layers.at(0)->GetYSizeByBlock(blockNum);
		int blockrow = EDS->Layers.at(0)->baseRef->getBlockRows();


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
			iX = Xmin + iCol*Cellsize + halfcell;
#ifdef BLOCKIO
			iY = Ymin + (ySize - blockNum*blockYsize - iRow - 1)*Cellsize + halfcell;
#else
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			iY = Ymin + (ySize - EDS->LayerRef->_pMetaData->_MBR.minIRow() - iRow - 1)*Cellsize + halfcell;
#endif

			bool isCal = true;	
			for (int k = 0; k < EDS->LayerSize; k++) {
				if (envVals[k] - EDS->Layers.at(k)->NoDataValue < VERY_SMALL)
					isCal = false;
			}
			if (!isCal) {
				data_prediction[i] = EDS->NoDataValue;//.push_back(EDS->NoDataValue);
				data_uncertainty[i] = EDS->NoDataValue;//.push_back(EDS->NoDataValue);
				continue;
			}

			vector<double> tmp;
			//envSimi[i] = new double[sampleCount];

			double* envSimi = new double[sampleCount];
			int tmpCount = 0;
			double simi_max = 0.0;
			for (int j = 0; j < sampleCount; j++)
			{
				se = SampleEnvUnits[j];
				sX = se->Loc->X;
				sY = se->Loc->Y;

				envSimi[j] = CalcSimi(se, envVals, EuclideanDistance);
				if (envSimi[j] > thred_envsimi) {
					double x_distance = sX - iX;
					double y_distance = sY - iY;
					tmp.push_back(sqrt(x_distance*x_distance + y_distance*y_distance)); //* 100000);//*100000);* 111*111*0.866
					tmpCount++;	
				}
				if (envSimi[j] > simi_max) {
					simi_max = envSimi[j];
				}
			}
						
			double sum1 = 0;	// sum of soil property * environmental similarity
			double sum2 = 0;	// sum of environmental similarities
			if (tmpCount > 0) {

				double maxdis, mindis;
				double adaR = AdaptiveR(&tmp, expR, &maxdis, &mindis);
				double m = Trimember(adaR, maxdis, mindis, 0.5, 1.0);
				tmpCount = 0;
				for (int j = 0; j < sampleCount; j++)
				{
					se = SampleEnvUnits[j];
					if (envSimi[j] > thred_envsimi) {

						simi_dis = pow(tmp[tmpCount], -m);
						sum1 += envSimi[j] * se->SoilVariable*simi_dis;
						sum2 += envSimi[j]*simi_dis;
						tmpCount++;
					}
				}
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