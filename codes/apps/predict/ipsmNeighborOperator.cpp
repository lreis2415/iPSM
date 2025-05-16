#include "ipsmNeighborOperator.h"
#include <iomanip>


// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

void solim::ipsmNeighborOperator::Initialize()
{
    this->thred_envsimi = 0.5;
	this->val_cannot_pred = -1.0; 
    this->nbrMin = 1;
    this->nbrMax = 100;
    this->Map_Prediction = nullptr;
	this->Map_Uncertainty = nullptr;
	this->n_step = 10;
	this->n_size = 10;
	this->n_step_amplifier = 1;
	this->n_range = 100;
}

bool solim::ipsmNeighborOperator::PredictMap_Property(double alpha){    
    this->n_step_amplifier = int(this->n_step / this->EDS->CellSize + 0.5);
    if (this->n_step_amplifier < 1) {
        this->n_step_amplifier = 1;
        this->n_size = int(this->n_range / this->EDS->CellSize + 0.5);
        this->n_size = this->n_size > 1 ? this->n_size : 1;
    }
    int count = this->EDS->TotalX * this->EDS->TotalY;
    int layerNum = this->EDS->Layers.size();
    if (nullptr == this->Map_Prediction) this->Map_Prediction = new EnvLayer();
    if (nullptr == this->Map_Uncertainty) this->Map_Uncertainty = new EnvLayer();
    this->Map_Prediction->CopyFrame(EDS->Layers[0]);
    this->Map_Uncertainty->CopyFrame(EDS->Layers[0]);

    float* data_prediction = this->Map_Prediction->EnvData;
    float* data_uncertainty = this->Map_Uncertainty->EnvData;

    int sampleCount = this->SampleEnvUnits.size();
    int segementCount = 1000;
    EnvUnit* se = nullptr;
    for (int i = 0; i < count; i++) {
        if (i % (count / segementCount) == 0)
        {
            //cout<<'\r';
            //cout<<"Completed "<<setw(5)<<(int((i*1.0/count+0.5/segementCount)*segementCount))/(segementCount/100.0)<<"%";
        }
        float* envVals = new float[EDS->LayerSize];
        EDS->GetEnvUnitValues(i, envVals);
        bool isCal = true;
        for (int k = 0; k < EDS->LayerSize; k++) {
            if (fabs(envVals[k] - EDS->Layers.at(k)->NoDataValue) < VERY_SMALL)
                isCal = false;
        }
        if (!isCal) {
            data_prediction[i] = this->EDS->NoDataValue;
            data_uncertainty[i] = -1;
            continue;
        }

        double simi_max = 0.0;
        double sum1 = 0;	// 求和（土壤属性值 * 环境相似度值）
        double sum2 = 0;	// 求和（环境相似度值）
        /// TODO: 重复定义了count. by lj.
        int count = 0;		// 计数，满足推测条件的样点数量
        for (int j = 0; j < sampleCount; j++)
        {
            se = this->SampleEnvUnits[j];
            if (!se->IsCal) {
                continue;
            }
            double envSimi = CalNbrSimi_Sample(i, envVals, se, alpha);//CalcSimi(se, envVals);

            if (envSimi == this->EDS->NoDataValue) {
                continue;
            }
            if (envSimi > simi_max) {
                simi_max = envSimi;
            }
            if (envSimi > this->thred_envsimi)			// 筛选出与待推测点环境相似度高的样点
            {
                sum1 += envSimi * se->SoilVariable;
                sum2 += envSimi;
                count++;
            }
        }

        //cout << i << ", " << simi_max << endl;

        if (count > 0)
        {
            data_prediction[i] = (1.0 * sum1 / sum2);
            //cout << i << ", " << simi_max << ", "<< 1.0 * sum1 / sum2 << endl;
            data_uncertainty[i] = (1 - simi_max);
        }
        else	// 没有可推测点，设为-1
        {
            data_prediction[i] = (this->val_cannot_pred);

            data_uncertainty[i] = (1 - simi_max);
        }
    }
    //cout<<'\r';
    //cout<<"Completed "<<setw(5)<<100.0<<"%";

    // 生成三个推测图层
    this->Map_Prediction->EnvData = data_prediction;
    this->Map_Uncertainty->EnvData = data_uncertainty;
    this->Map_Prediction->NoDataValue = this->EDS->NoDataValue;
    this->Map_Uncertainty->NoDataValue = -1;
    this->Map_Prediction->CalcStat();
    this->Map_Uncertainty->CalcStat();
    this->Map_Prediction->Writeout(map_predict_filename);
    this->Map_Uncertainty->Writeout(map_uncer_filename);
    se = nullptr;
    return true;
}

bool solim::ipsmNeighborOperator::validate(double alpha,vector<EnvUnit*>validPoints) {
	int layerNum = this->EDS->Layers.size();
	int validPoint_num = validPoints.size();
	float *data_prediction = new float[validPoint_num];
	float *data_uncertainty = new float[validPoint_num];
	int count_valid = 0;
	double err_sum = 0;
	double err_sqr_sum = 0;

	int sampleCount = this->SampleEnvUnits.size();
	int segementCount = 1;
	EnvUnit *se = nullptr;
	for (int i = 0; i < validPoints.size(); i++) {
		if (i % (validPoint_num / segementCount) == 0)
		{
			cout << '\r';
			cout << "Completed " << setw(5) << (int((i*1.0 / validPoint_num + 0.5 / segementCount)*segementCount)) / (segementCount / 100.0) << "%";
		}
		float *envVals = new float[EDS->LayerSize];
		for (int f = 0; f < EDS->LayerSize; f++) {
			envVals[f] = validPoints[i]->EnvValues[f];
		}
		long ind = validPoints[i]->Loc->Row*EDS->TotalX + validPoints[i]->Loc->Col;
		bool isCal = true;
		for (int k = 0; k < EDS->LayerSize; k++) {
			if (fabs(envVals[k] - EDS->Layers.at(k)->NoDataValue) < VERY_SMALL)
				isCal = false;
		}
		if (!isCal) {
			data_prediction[i] = this->EDS->NoDataValue;
			data_uncertainty[i] = -1;
			continue;
		}

		double simi_max = 0.0;
		double sum1 = 0;	// 求和（土壤属性值 * 环境相似度值）
		double sum2 = 0;	// 求和（环境相似度值）
							/// TODO: 重复定义了count. by lj.
		int count = 0;		// 计数，满足推测条件的样点数量
		for (int j = 0; j < sampleCount; j++)
		{
			se = this->SampleEnvUnits[j];
			if (!se->IsCal) {
				continue;
			}
			double envSimi = CalNbrSimi_Sample(ind, envVals, se, alpha);

			if (envSimi == this->EDS->NoDataValue) {
				continue;
			}
			if (envSimi > simi_max) {
				simi_max = envSimi;
			}
			if (envSimi > this->thred_envsimi)			// 筛选出与待推测点环境相似度高的样点
			{
				sum1 += envSimi * se->SoilVariable;
				sum2 += envSimi;
				count++;
			}
		}

		//cout << i << ", " << simi_max << endl;

		if (count > 0)
		{
			data_prediction[i] = (1.0 * sum1 / sum2);
			double err = validPoints[i]->SoilVariable - data_prediction[i];
			err_sum += fabs(err);
			err_sqr_sum += err*err;
			count_valid++;
			//cout << i << ", " << simi_max << ", "<< 1.0 * sum1 / sum2 << endl;
			data_uncertainty[i] = (1 - simi_max);
		}
		else	// 没有可推测点，设为-1
		{
			data_prediction[i] = (this->val_cannot_pred);

			data_uncertainty[i] = (1 - simi_max);
		}
	}
	cout << '\r';
	cout << "Completed " << setw(5) << 100.0 << "%";
	cout << '\n' << "RMSE: " << sqrt(err_sqr_sum / count_valid);
	cout << '\n' << "MAE: " << sqrt(err_sum / count_valid);
	se = nullptr;
	return true;
}


double solim::ipsmNeighborOperator::CalNbrSimi_Sample(long pixelIndex, float *pixel,
    EnvUnit *sample, double alpha){   
    int sampleRow = sample->Loc->Row;
    int sampleCol = sample->Loc->Col;    
    double nbrSimiSample = 1;
    for(int f = 0; f < EDS->LayerSize; ++f){
        double nbrSimiVariable = CalNbrSimi_Variable(pixelIndex,
            sampleRow, sampleCol, f, pixel, sample, alpha);
        if(nbrSimiVariable == this->EDS->NoDataValue){
            return this->EDS->NoDataValue;
        }

        if(nbrSimiVariable < nbrSimiSample){
            nbrSimiSample = nbrSimiVariable;
        }
    }
    return nbrSimiSample;    
}	

double solim::ipsmNeighborOperator::CalNbrSimi_Variable(int indexPixel, 
    int sampleRow, int sampleCol, int f, float *pixel, EnvUnit *sample,
    double alpha){
	double lyr_mean = this->EDS->Layers[f]->Data_Mean;
	double lyr_std = this->EDS->Layers[f]->Data_StdDev;
	DataTypeEnum lyr_type = this->EDS->Layers[f]->DataType;
	long xSize = this->EDS->TotalX;
	long ySize = this->EDS->TotalY;
	double envPixel = pixel[f];
	double envSample = sample->EnvValues[f];
	//待推测点和样点基于点的相似性
	double pointSimilarity = CalcSimi_Single_Gaussian(envSample, envPixel, lyr_mean, lyr_std, lyr_type);
	//cout << pointSimilarity <<endl;

	//求中心栅格相似性和不同size的环形的相似性的权重
	double dist_w_sum = 1.0 / pow(1, alpha);

	//求邻域相似性
	double nbrSimilarity = pointSimilarity * dist_w_sum;
	int pCol = indexPixel % EDS->TotalX;
	int pRow = indexPixel / EDS->TotalX;
	for (int size = 1; size <= n_size; size++) {

		//待推测点和样点环形的环境条件向量
		vector<float> envVectPixel, envVectSample, envVectSampleMoved;
		vector<int> rCord, cCord;
		getAnnulusCord(size, rCord, cCord);
		this->getEnvVector(pCol, pRow, this->EDS->Layers[f], rCord, cCord, envVectPixel);
		this->getEnvVector(sample->Loc->Col, sample->Loc->Row, this->EDS->Layers[f], rCord, cCord, envVectSample);
		vector<int>().swap(rCord);
		vector<int>().swap(cCord);

		//求环形相似性
		double annulusSimi = CalcVectorSimi(envVectSample, envVectPixel, this->EDS->Layers[f]->NoDataValue, f);

		double w = 1.0 / pow(size * n_step + 1, alpha);
		nbrSimilarity += w * annulusSimi;
		dist_w_sum += w;

	}
	//求样点和待推测点单点相似性    
	//cout << nbrSimilarity << endl;
	return nbrSimilarity / dist_w_sum;
}



void solim::ipsmNeighborOperator::getEnvVector(int pCol, int pRow, EnvLayer *Layer, vector<int>rCord, 
    vector<int>cCord, vector<float> & envVector){
	for (int l = 0; l < rCord.size(); ++l) {
		int rowNbrPoint = pRow + rCord[l];
		int colNbrPoint = pCol + cCord[l];
		if (rowNbrPoint < 0 || rowNbrPoint >= Layer->YSize || colNbrPoint < 0 || colNbrPoint >= Layer->XSize) {
			continue;
		}
		float envNbr = Layer->EnvData[rowNbrPoint * Layer->XSize + colNbrPoint];
		envVector.push_back(envNbr);
	}

}

//calculate the similarity of two env vectors
double solim::ipsmNeighborOperator::CalcVectorSimi(vector<float> e1, vector<float> e2, double NoDataValue,int f) {
	double lyr_mean = this->EDS->Layers[f]->Data_Mean;
	double lyr_std = this->EDS->Layers[f]->Data_StdDev;
	DataTypeEnum lyr_type = this->EDS->Layers[f]->DataType;
	if (lyr_type == CATEGORICAL) {
		int mode1, mode2, unqNum1, unqNum2;
		getModeUnique(e1, mode1, unqNum1);
		getModeUnique(e2, mode2, unqNum1);
		if (mode1 == mode2) {
			int type_num = 14;
			if (f > 1) type_num = 6;
			return 1 - fabs(unqNum1 - unqNum1) / type_num;
		}
		else return 0;
	}
	else {
		int vectLength = e1.size();
		double mean1 = std::accumulate(e1.begin(), e1.end(), 0.0) / e1.size();
		double max1 = *max_element(e1.begin(), e1.end());
		double min1 = *min_element(e1.begin(), e1.end());
		double range1 = max1 - min1;
		double mean2 = std::accumulate(e2.begin(), e2.end(), 0.0) / e2.size();
		double max2 = *max_element(e2.begin(), e2.end());
		double min2 = *min_element(e2.begin(), e2.end());
		double range2 = max2 - min2;
		double mean_simi = CalcSimi_Single_Gaussian(mean1, mean2, lyr_mean, lyr_std, lyr_type);
		//double max_simi = CalcSimi_Single_Gaussian(max1, max2, lyr_mean, lyr_std, lyr_type);
		//double min_simi = CalcSimi_Single_Gaussian(min1, min2, lyr_mean, lyr_std, lyr_type);
		//return (mean_simi + max_simi + min_simi) / 3;
		double range_simi = 1 - fabs(range1 - range2) / 3.0 / lyr_std;
		range_simi = (range_simi < 0) ? 0 : range_simi;
		return 0.5 * mean_simi + 0.5 * range_simi;
	}
}

//get the relative coordinate of annulus pixels to the interest pixel
//when the outer radius of the annulus is size and the inner radius is size-1
//the default shape of the neighborhood is circle
void solim::ipsmNeighborOperator::getAnnulusCord(int size, vector<int>& rCord, vector<int>& cCord) {
    int outer2 = pow(n_step_amplifier * size, 2);             //外环半径平方
    int inner2 = pow(n_step_amplifier * (size - 1), 2); //内环半径平方
    for (int r = -n_step_amplifier * size - 1; r <= n_step_amplifier * size + 1; r++) {
        for (int c = -n_step_amplifier * size - 1; c <= n_step_amplifier * size + 1; c++) {
            if (r * r + c * c <= outer2 && r * r + c * c > inner2) {
                rCord.push_back(r);
                cCord.push_back(c);
            }
        }
    }
}

void solim::ipsmNeighborOperator::getModeUnique(vector<float> v, int& mode, int& unique_num) {
    vector<int> i_v;
    for (int k = 0; k < v.size(); k++) {
        i_v.push_back(int(v[k]));
    }
    sort(i_v.begin(), i_v.end());
    mode = i_v[0];
    int local_freq = 0;
    int final_freq = 0;
    unique_num = 0;
    for (vector<int>::iterator it = i_v.begin(); it != i_v.end() - 1; it++) {
        if (*it == *(it + 1)) local_freq++;
        else {
            local_freq = 0;
            unique_num++;
        }
        if (local_freq > final_freq)
        {
            final_freq = local_freq;
            mode = *it;
        }
    }
}
