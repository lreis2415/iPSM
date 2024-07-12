/*!
 * @brief IPSMFILTER algorithm.
 * @details 对于算法，在头文件文件开头添加标准格式信息，包含但不限于以下几点：by fnq, 20210905.
 *          1. 算法原理简介
 *          2. 参考文献
 *          3. 算法作者及开发者
 * @todo 添加格式化代码注释. by fnq..
 * @version 1.0
 * @revision  21-09-05 fnq - initial version
 */
#ifndef IPSMFILTEROPERATOR_H_
#define IPSMFILTEROPERATOR_H_

#include "CommonOperator.h"
#include <iostream>
#include <iomanip>

using namespace std;

namespace solim
{
	class ipsmFilterOperator : public CommonOperator
	{
		// constructor functions
	public:
		ipsmFilterOperator() : CommonOperator()
		{
			this->Initialize();
		}

        // TODO, use explicit keyword for singe-parameter function. by lj.
		explicit ipsmFilterOperator(EnvDataset *envDataset) : CommonOperator(envDataset)
		{
			this->Initialize();
		}

		ipsmFilterOperator(EnvDataset *envDataset, vector<EnvUnit *> sampleEnvUnits) :
            CommonOperator(envDataset, sampleEnvUnits)
		{
			this->Initialize();
		}

		~ipsmFilterOperator(void) {}; 
		// methods
	public:
		void Initialize();
		bool PredictMap_Property();
		double CalcSimi(EnvUnit *se, EnvUnit *e, SimilarityTypeEnum simitype = GaussianDistance);
		double CalcSimi(EnvUnit *se, float *envVals, SimilarityTypeEnum simitype = GaussianDistance);
		double CalcUncertainty(EnvUnit* e, SimilarityTypeEnum simitype = GaussianDistance);
		double CalcUncertainty(float* envVals, SimilarityTypeEnum simitype = GaussianDistance);
		double CalcFilterUncer(float* envVals);
		double CalcUncertainty(double simi_max, float *envVals);
		void setDataFactors(vector<string> &dataFactorList);

		// fields
	public:
		double thred_envsimi;
		double val_cannot_pred;
		EnvLayer *Map_Prediction;
		EnvLayer *Map_Uncertainty;  
		string map_predict_filename;
		string map_uncer_filename;
		int Count_Geo, Count_Terrain, Count_Climate, Count_Vege, Count_Other;
		vector<DataFactorEnum> dataFactors;
	};
}
#endif
