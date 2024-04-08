#include "ipsmFilterOperator.h"


using namespace std;

void solim::ipsmFilterOperator::Initialize()
{
	this->thred_envsimi = 0.5;
	this->val_cannot_pred = -1.0;
	this->Map_Prediction = new EnvLayer();
	this->Map_Uncertainty = new EnvLayer();
}

bool solim::ipsmFilterOperator::PredictMap_Property()
{
	int xSize = this->EDS->TotalX;
	int ySize = this->EDS->TotalY;
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
			se = this->SampleEnvUnits[0];
			int null_num = 0;
			for (int k = 0; k < EDS->LayerSize; k++) {
				if (envVals[k] - EDS->Layers.at(k)->NoDataValue < VERY_SMALL)
					null_num++;
			}
			if (null_num == se->EnvValues.size())
			{
				data_prediction[i] = EDS->NoDataValue;
				data_uncertainty[i] = EDS->NoDataValue;
				continue;
			}
			double simi_max = 0.0;
			double sum1 = 0;	// 求和（土壤属性值 * 环境相似度值）
			double sum2 = 0;	// 求和（环境相似度值）

			int count = 0;		// 计数，满足推测条件的样点数量
			for (int j = 0; j < sampleCount; j++)
			{
				se = this->SampleEnvUnits[j];
				double envSimi = CalcSimi(se, envVals);//此处采用了FilteriPSM的相似度计算函数
				if (envSimi > simi_max)
				{
					simi_max = envSimi;
				}if (envSimi > this->thred_envsimi)			// 筛选出与待推测点环境相似度高的样点
				{
					sum1 += envSimi * se->SoilVariable;
					sum2 += envSimi;
					count++;
				}
			}
			double uncer = CalcUncertainty(simi_max, envVals);
			if (count > 0)
			{
				data_prediction[i] = 1.0 * sum1 / sum2;
				data_uncertainty[i] = uncer;
			}
			else	// 没有可推测点，设为-1
			{
				data_prediction[i] = EDS->NoDataValue;
				data_uncertainty[i] = uncer;
			}	
		}
		cout<<'\r';
		cout<<"Completed "<<setw(5)<<100.0<<"%";

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

/*!
* @brief 计算两个环境单元之间的综合相似度(也可以计算当位置某个环境变量为空值时的相似度)
* @param se
* @param e
* @return
* by fnq 20210912
* changed to override function by zfh 2021/12/1
*/
double solim::ipsmFilterOperator::CalcSimi(EnvUnit *se, EnvUnit *e, SimilarityTypeEnum simitype) {
	int envSize = e->EnvValues.size();
	float* envVals = new float[envSize];
	for (int i = 0; i < envSize; i++) {
		envVals[i] = e->EnvValues[i];
	}
	return CalcSimi(se, envVals, simitype);
}

double solim::ipsmFilterOperator::CalcSimi(EnvUnit *se, float *envVals, SimilarityTypeEnum simitype){
	double simi = -1.;
	int null_num = 0;
	for (int i = 0; i < se->EnvValues.size(); i++) {
		float nodata = EDS->Layers.at(i)->NoDataValue;
		if (fabs(envVals[i] - nodata) < VERY_SMALL) {
			null_num++;
		}
	}
	if (null_num == se->EnvValues.size()) {
		return -1.;    // 不参与计算的点
	}
	if (se->EnvValues.size() == EDS->LayerSize) {
		if (simitype == GaussianDistance) {
			simi = 1;
			for (int i = 0; i < se->EnvValues.size(); i++) {
				double mean = EDS->Layers.at(i)->Data_Mean;
				double stdDev = EDS->Layers.at(i)->Data_StdDev;
				DataTypeEnum dataType = se->DataTypes.at(i);
				double simi_temp = CalcSimi_Single_Gaussian(se->EnvValues.at(i), envVals[i], mean, stdDev, dataType);
				if (simi_temp < 0.0) {
					simi_temp = 1.0;
				}
				if (simi_temp < simi) {
					simi = simi_temp;
				}
			}
			return simi;
		}
		else if (simitype == EuclideanDistance) {
			simi = 1;
			for (int i = 0; i < se->EnvValues.size(); i++) {
				double range = EDS->Layers.at(i)->Data_Range;
				DataTypeEnum dataType = se->DataTypes.at(i);
				double simi_temp = CalcSimi_Single_Euclidean(se->EnvValues.at(i), envVals[i], range, dataType);
				if (simi_temp < 0.0) {
					simi_temp = 1.0;
				}
				if (simi_temp < simi) {
					simi = simi_temp;
				}
			}
			return simi;
		}
	}
	else {
		return -1.;
	}
}

void solim::ipsmFilterOperator::setDataFactors(vector<string> &dataFactorList){
	Count_Geo = 0;
	Count_Terrain = 0;
	Count_Climate = 0;
	Count_Vege = 0;
	Count_Other = 0;
	dataFactors.clear();
	dataFactors.shrink_to_fit();
	for (int i = 0; i < dataFactorList.size(); i++) {
		string datafactor = dataFactorList[i];
		transform(datafactor.begin(), datafactor.end(), datafactor.begin(), ::toupper);
		if (datafactor == "GEOLOGY") {
			dataFactors.push_back(GEOLOGY);
			Count_Geo = Count_Geo + 1;
		}
		else if (datafactor == "TERRAIN") {
			dataFactors.push_back(TERRAIN);
			Count_Terrain = Count_Terrain + 1;
		}
		else if (datafactor == "CLIMATE") {
			dataFactors.push_back(CLIMATE);
			Count_Climate = Count_Climate + 1;
		}
		else if (datafactor == "VEGETATION") {
			dataFactors.push_back(VEGETATION);
			Count_Vege = Count_Vege + 1;
		}
		else {
			dataFactors.push_back(OTHERS);
			Count_Other = Count_Other + 1;
		}
	}
}


double solim::ipsmFilterOperator::CalcUncertainty(EnvUnit *e, SimilarityTypeEnum simitype) {
	int envSize = e->EnvValues.size();
	float* envVals = new float[envSize];
	for (int i = 0; i < envSize; i++) {
		envVals[i] = e->EnvValues[i];
	}
	return CalcUncertainty(envVals, simitype);
}

double solim::ipsmFilterOperator::CalcUncertainty(float *envVals, SimilarityTypeEnum simitype) {
	double simi = 0;
	double simi_temp = 0;
	for (vector<EnvUnit *>::iterator it = SampleEnvUnits.begin(); it != SampleEnvUnits.end(); ++it) {
		simi_temp = CalcSimi(*it, envVals, simitype);
		if (simi < simi_temp) {
			simi = simi_temp;
		}
	}
	double filterUncer = CalcFilterUncer(envVals);
	return 1 - simi_temp + filterUncer - (1 - simi_temp)*filterUncer;
}

double solim::ipsmFilterOperator::CalcFilterUncer(float* envVals) {
	double ipsmFilterUncer_temp = 0.;
	//double simi = 0.;
	//double simi_temp = 0.;
	double uncertainty = 0.;
	int count_factor = 0; //统计环境变量种类数，即输入的环境变量属于五大成土类型中的几类，如输入的是高程、坡度、降水，此时即为2种类型（地形类、气候类）

						  //统计每个类型的环境变量中有几个环境变量空值
	int count_filtergeo = 0, count_filterterrain = 0, count_filterclimate = 0, count_filtervege = 0, count_filterother = 0;
	int nodatavalue = this->EDS->NoDataValue;

	// 统计环境变量种类数，
	if (Count_Geo != 0) count_factor = count_factor + 1;
	if (Count_Terrain != 0) count_factor = count_factor + 1;
	if (Count_Climate != 0) count_factor = count_factor + 1;
	if (Count_Vege != 0) count_factor = count_factor + 1;
	if (Count_Other != 0) count_factor = count_factor + 1;

	for (int i = 0; i < dataFactors.size() && i < EDS->LayerSize; i++) {
		float nodata = EDS->Layers.at(i)->NoDataValue;
		if (fabs(envVals[i] - nodata) < VERY_SMALL) {
			//统计各类中的空值环境变量的个数
			DataFactorEnum dataFactor = dataFactors.at(i);
			if (dataFactor == GEOLOGY) count_filtergeo = count_filtergeo + 1;
			else if (dataFactor == TERRAIN) count_filterterrain = count_filterterrain + 1;
			else if (dataFactor == CLIMATE) count_filterclimate = count_filterclimate + 1;
			else if (dataFactor == VEGETATION) count_filtervege = count_filtervege + 1;
			else count_filterother = count_filterother + 1;
		}
		else {
			ipsmFilterUncer_temp = 0.;
		}
	}
	if(dataFactors.size()==0) return 1-ipsmFilterUncer_temp;

	//计算过滤空值环境变量后的推测不确定性
	if (Count_Geo != 0) ipsmFilterUncer_temp = ipsmFilterUncer_temp + (double)count_filtergeo / (count_factor * Count_Geo);
	if (Count_Terrain != 0) ipsmFilterUncer_temp = ipsmFilterUncer_temp + (double)count_filterterrain / (count_factor * Count_Terrain);
	if (Count_Climate != 0) ipsmFilterUncer_temp = ipsmFilterUncer_temp + (double)count_filterclimate / (count_factor * Count_Climate);
	if (Count_Vege != 0) ipsmFilterUncer_temp = ipsmFilterUncer_temp + (double)count_filtervege / (count_factor * Count_Vege);
	if (Count_Other != 0) ipsmFilterUncer_temp = ipsmFilterUncer_temp + (double)count_filterother / (count_factor * Count_Other);

	return ipsmFilterUncer_temp;
}

double solim::ipsmFilterOperator::CalcUncertainty(double simi_max, float *envVals) {
	double filterUncer = CalcFilterUncer(envVals);
	return 1 - simi_max + filterUncer - (1 - simi_max)*filterUncer;
}