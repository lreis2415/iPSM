#include "ipsmNeighborOperator.h"
#include "solim_math.h"
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
}


bool solim::ipsmNeighborOperator::CalcChrScale(vector<string> map_scale_names){
	this->Map_scale = new EnvDataset;

    int layNum = this->EDS->Layers.size();
    int xSize = this->EDS->TotalX;
    int ySize = this->EDS->TotalY;
    int count = xSize * ySize;


    for(int f = 0; f < layNum; ++f){//单个变量
        this->EDS->Layers[f]->CalcStat();
        EnvLayer *scaleLayer = new EnvLayer();        
		scaleLayer->CopyFrame(EDS->Layers[0]);
        float* EnvScales = scaleLayer->EnvData;
        
		double lyr_mean = this->EDS->Layers[f]->Data_Mean;
		double lyr_std = this->EDS->Layers[f]->Data_StdDev;
		DataTypeEnum lyr_type = this->EDS->Layers[f]->DataType;
              
        for (int i = 0; i < count; i++){
             double envValue = this->EDS->Layers[f]->EnvData[i];
             int row =i/xSize;
             int col = i%xSize;        
                     
            //逐栅格计算特征邻域 
             if(envValue == this->EDS->Layers[f]->NoDataValue){
                EnvScales[i]=(this->EDS->NoDataValue);
                continue;
            }

            int size = nbrMin, flag = 1;
            int chrSize = 0, annulusIndex = 0, diffIndex = 0;
            vector <double> annulusSimilarity, similarityDiff;
            double sumPointSimilarity = 0;
            while(flag == 1 && size < nbrMax){
                vector<int> rCord, cCord;
                getAnnulusCord(size, rCord, cCord);//顺时针记录环形内的栅格与中心栅格的相对坐标
                sumPointSimilarity = 0;
                int cnt = 0;            
                for(int n = 0; n < rCord.size(); ++n){
                    int rowNeighbor = row + rCord[n];
                    if(rowNeighbor < 0 || rowNeighbor >= ySize)
                        continue;
                    int colNeighbor = col + cCord[n];
                    if(colNeighbor < 0 || colNeighbor >= xSize)
                        continue;
                    double envNeighbor = this->EDS->Layers[f]->EnvData[rowNeighbor * xSize + colNeighbor];
                    if(envNeighbor == this->EDS->Layers[f]->NoDataValue){
                        continue; 
                    }
                                       
					double pointSimilarity = CalcSimi_Single_Gaussian(envNeighbor, envValue, lyr_mean, lyr_std, lyr_type);
                    sumPointSimilarity += pointSimilarity;	
                    cnt++;						
                }
                rCord.clear();
                cCord.clear();
                if(cnt == 0){
                    flag = 0;
                    break;
                }
                //cout << size << ", " << sumPointSimilarity << ", " << cnt << endl;
                annulusSimilarity.push_back(sumPointSimilarity / cnt);		
                if(annulusIndex == 0){			
                    annulusIndex++;	
                    size++;	
                    continue;
                }			
                diffIndex = annulusIndex - 1;
                similarityDiff.push_back(annulusSimilarity[annulusIndex - 1] - annulusSimilarity[annulusIndex]); 				
                if(diffIndex < 2){
                    annulusIndex++;	
                    size++;	
                    continue;
                }
                if(similarityDiff[diffIndex - 1] > similarityDiff[diffIndex - 2] && similarityDiff[diffIndex - 1] > similarityDiff[diffIndex]){
                    chrSize = diffIndex;
                    flag = 0;								
                }
                annulusIndex++;
                size++;	
                vector<int>().swap(rCord);
                vector<int>().swap(cCord);            
            }

            //nbrMin和nbrMax之间没找到拐点的情况
            if(size == nbrMax){
                double envMean = 0, envSum = 0, varSum = 0, difMean = 0;
                int cnt = 0;
                for(int r = -nbrMax; r <= nbrMax; ++r){
                    int rowNeighbor = row + r;
                    if(rowNeighbor < 0 || rowNeighbor >= ySize)
                        continue;
                    for(int c = -nbrMax; c <= nbrMax; ++c){
                        int colNeighbor = col + c;
                        if(colNeighbor < 0 || colNeighbor >= xSize)
                            continue;
                        double envNeighbor = this->EDS->Layers[f]->EnvData[rowNeighbor * xSize + colNeighbor];
                        if(envNeighbor == this->EDS->Layers[f]->NoDataValue)
                            continue;
                        envSum += envNeighbor;
                        cnt++;						
                    }				
                }
                envMean = envSum / cnt;
                for(int r = -nbrMax; r <= nbrMax; ++r){
                    int rowNeighbor = row + r;
                    if(rowNeighbor < 0 || rowNeighbor >= ySize)
                        continue;
                    for(int c = -nbrMax; c <= nbrMax; c++){
                        int colNeighbor = col + c;
                        if(colNeighbor < 0 || colNeighbor >= xSize)
                            continue;
                        double envNeighbor = this->EDS->Layers[f]->EnvData[rowNeighbor * xSize + colNeighbor];
                        if(envNeighbor == this->EDS->Layers[f]->NoDataValue)
                            continue;
                        varSum += pow(envNeighbor - envMean, 2);						
                    }				
                }
                difMean = varSum / cnt;		
                if(difMean >= this->EDS->Layers[f]->Data_StdDev){
                    chrSize = nbrMin;
                }else{
                    chrSize = nbrMax;
                }		
            }
            vector<double>().swap(annulusSimilarity);
            vector<double>().swap(similarityDiff);
            //cout << chrSize << endl;
            EnvScales[i]=(chrSize);                      
        }

        scaleLayer->EnvData = EnvScales;
        scaleLayer->NoDataValue = this->EDS->NoDataValue;
        this->Map_scale->AddLayer(scaleLayer);
		if (map_scale_names.size() == EDS->LayerSize)
			scaleLayer->Writeout(map_scale_names[f]);
    }
    return true;
}

bool solim::ipsmNeighborOperator::PredictMap_Property(double alpha){    
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
    EnvUnit *se = nullptr;    
    for (int i = 0; i < count; i++){
        if(i % (count/segementCount) == 0)
        {
            cout<<'\r';
            cout<<"Completed "<<setw(5)<<(int((i*1.0/count+0.5/segementCount)*segementCount))/(segementCount/100.0)<<"%";
        }
		float *envVals = new float[EDS->LayerSize];
		EDS->GetEnvUnitValues(i, envVals);
		bool isCal = true;
		for (int k = 0; k < EDS->LayerSize; k++) {
			if (fabs(envVals[k] - EDS->Layers.at(k)->NoDataValue) < VERY_SMALL)
				isCal = false;
		}
        if(!isCal){
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
            if(!se->IsCal){
                continue;
            }
            double envSimi = CalNbrSimi_Sample(i, envVals, se, alpha);

            if(envSimi == this->EDS->NoDataValue){
                continue;
            }
            if (envSimi > simi_max){
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
        } else	// 没有可推测点，设为-1
        {
            data_prediction[i] = ( this->val_cannot_pred );
            
            data_uncertainty[i] = (1 - simi_max);
        }	
    }
    cout<<'\r';
    cout<<"Completed "<<setw(5)<<100.0<<"%";

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
        //综合特征邻域
        int pixelScale = this->Map_scale->Layers[f]->EnvData[indexPixel];
        int sampleIndex = sampleRow * xSize + sampleCol;
        int sampleScale = this->Map_scale->Layers[f]->EnvData[sampleIndex];
        double envPixel = pixel[f];
        double envSample = sample->EnvValues[f];
        if(pixelScale == this->EDS->NoDataValue || sampleScale == this->EDS->NoDataValue){
            return this->EDS->NoDataValue;
        }
        int scaleIntegration = (pixelScale<sampleScale)?pixelScale:sampleScale;
        //待推测点和样点基于点的相似性
		double pointSimilarity = CalcSimi_Single_Gaussian(envSample, envPixel, lyr_mean, lyr_std, lyr_type);
        //cout << pointSimilarity <<endl;

        //求中心栅格相似性和不同size的环形的相似性的权重
        double *weight = new double[scaleIntegration + 1];
        CalcAnnulusWeight(weight, scaleIntegration, alpha); 

        //求邻域相似性
        double nbrSimilarity = 0;      
		int pCol = indexPixel%EDS->TotalX;
		int pRow = indexPixel / EDS->TotalX;
        for(int size = 1; size <= scaleIntegration; size++){

            //待推测点和样点环形的环境条件向量
            vector<float> envVectPixel, envVectSample, envVectSampleMoved;
            vector<int> rCord, cCord;
            getAnnulusCord(size, rCord, cCord);
            this->getEnvVector(pCol,pRow, this->EDS->Layers[f], rCord, cCord, envVectPixel);
            this->getEnvVector(sample->Loc->Col, sample->Loc->Row, this->EDS->Layers[f], rCord, cCord, envVectSample);
            vector<int>().swap(rCord);
            vector<int>().swap(cCord);  

            //求环形相似性
            double annulusSimi = CalcVectorSimi(envVectPixel, envVectSample, this->EDS->Layers[f]->NoDataValue);
            for(int n = 1; n < envVectSample.size(); ++n){            
                for(int m = 0; m < envVectSample.size() - n; ++m){
                    envVectSampleMoved.push_back(envVectSample[m + n]);

                }
                for(int m = envVectSample.size() - n; m < envVectSample.size(); ++m){
                    envVectSampleMoved.push_back(envVectSample[m - (envVectSample.size() -n)]);
                }
                double vectSimiMoved = CalcVectorSimi(envVectPixel, envVectSampleMoved, 
                    this->EDS->Layers[f]->NoDataValue);
                if(annulusSimi < vectSimiMoved){
                    annulusSimi = vectSimiMoved;
                }            
                envVectSampleMoved.clear();
            }
            //cout << annulusSimi << ", ";
            vector<float>().swap(envVectPixel);
            vector<float>().swap(envVectSample);
            vector<float>().swap(envVectSampleMoved);
            nbrSimilarity += weight[size] * annulusSimi;
        }
        //求样点和待推测点单点相似性    

        nbrSimilarity += weight[0] * pointSimilarity;
        //cout << nbrSimilarity << endl;
        delete [] weight;
        return nbrSimilarity;
}



void solim::ipsmNeighborOperator::getEnvVector(int pCol, int pRow, EnvLayer *Layer, vector<int>rCord, 
    vector<int>cCord, vector<float> & envVector){
        for(int l = 0; l < rCord.size(); ++l){
            int rowNbrPoint = pRow + rCord[l];
            int colNbrPoint = pCol + cCord[l];
            if(rowNbrPoint < 0 || rowNbrPoint >= Layer->YSize||colNbrPoint < 0 || colNbrPoint >= Layer->XSize){
                envVector.push_back(Layer->NoDataValue);
                continue;
            }
            float envNbr = Layer->EnvData[rowNbrPoint * Layer->XSize + colNbrPoint];
            envVector.push_back(envNbr);
        }

}

//calculate the similarity of two env vectors
double solim::ipsmNeighborOperator::CalcVectorSimi(vector<float> e1, vector<float> e2, double NoDataValue) {
	double innerProduct = 0;
	double modeSquare1 = 0;
	double modeSquare2 = 0;
	int vectLength = e1.size();
	for (int l = 0; l < vectLength; ++l) {
		if (e1[l] == NoDataValue || e2[l] == NoDataValue)
			continue;
		innerProduct += e1[l] * e2[l];
		modeSquare1 += e1[l] * e1[l];
		modeSquare2 += e2[l] * e2[l];
	}
	if (innerProduct == 0)
		return 0;
	double simi = innerProduct / sqrt(modeSquare1) / sqrt(modeSquare2);
	return simi;
}

//get the relative coordinate of annulus pixels to the interest pixel
//when the outer radius of the annulus is size and the inner radius is size-1
//the default shape of the neighborhood is circle
void solim::ipsmNeighborOperator::getAnnulusCord(int size, vector<int>& rCord, vector<int>& cCord) {
	int outer2 = size * size;             //外环半径平方
	int inner2 = (size - 1) * (size - 1); //内环半径平方
	vector<int> rDist;
	vector<int> cDist;
	vector<int> rDistDescent; //rDist从小到大排序
	vector<int> cDistDescent;
	for (int r = size; r >= 0; r--) {
		for (int c = size - r; c <= size; c++) {
			if (r * r + c * c <= outer2 && r * r + c * c > inner2) {
				rDist.push_back(r);
				cDist.push_back(c);
			}
		}
	}
	for (int v = 0; v < rDist.size(); v++) {
		int r = rDist[rDist.size() - v - 1];
		int c = cDist[cDist.size() - v - 1];
		rDistDescent.push_back(r);
		cDistDescent.push_back(c);
	}
	for (int v = 0; v < rDist.size(); v++) {
		int r = -rDist[v];
		int c = cDist[v];
		rCord.push_back(r);
		cCord.push_back(c);
	}
	for (int v = 0; v < rDist.size(); v++) {
		int r = rDistDescent[v];
		int c = cDistDescent[v];
		if (r != rCord[rCord.size() - 1] || c != cCord[cCord.size() - 1]) {
			rCord.push_back(r);
			cCord.push_back(c);
		}
	}
	for (int v = 0; v < rDist.size(); v++) {
		int r = rDist[v];
		int c = -cDist[v];
		if (r != rCord[rCord.size() - 1] || c != cCord[cCord.size() - 1]) {
			rCord.push_back(r);
			cCord.push_back(c);
		}
	}
	for (int v = 0; v < rDist.size(); v++) {
		int r = -rDistDescent[v];
		int c = -cDistDescent[v];
		if ((r != rCord[rCord.size() - 1] || c != cCord[cCord.size() - 1]) && (r != rCord[0] || c != cCord[0])) {
			rCord.push_back(r);
			cCord.push_back(c);
		}
	}
	vector<int>().swap(rDist);
	vector<int>().swap(cDist);
	vector<int>().swap(rDistDescent);
	vector<int>().swap(cDistDescent);
}

void solim::ipsmNeighborOperator::CalcAnnulusWeight(double* weight, int scale, double alpha) {
	int length = scale + 1;
	double* distance = new double[length];
	double sum = 0;
	for (int d = 0; d < length; d++) {
		distance[d] = pow(d + 0.5, alpha);
		sum += 1 / distance[d];
	}
	for (int d = 0; d < length; d++) {
		weight[d] = 1 / distance[d] / sum;
		//cout << weight[d] << endl;
	}
	delete[] distance;
}