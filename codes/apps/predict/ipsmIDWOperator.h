/*!
 * @brief IPSM algorithm.
 * @details 对于算法，在头文件文件开头添加标准格式信息，包含但不限于以下几点：by lj.
 *          1. 算法原理简介
 *          2. 参考文献
 *          3. 算法作者及开发者
 * @todo
 * @version 1.0
 * @revision  17-11-21 zhanglei - initial version
 *            17-11-27 lj       - code revision and format
 */
#ifndef IPSMIDWOPERATOR_H_
#define IPSMIDWOPERATOR_H_

#include "CommonOperator.h"

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

namespace solim {
class ipsmIDWOperator: public CommonOperator {
public:
    ipsmIDWOperator():
		thred_envsimi(0.5), val_cannot_pred(NODATA), r(0.5),
        Map_Prediction(nullptr), Map_Uncertainty(nullptr) {
    }

    explicit ipsmIDWOperator(EnvDataset* envDataset)
        : CommonOperator(envDataset), 
		thred_envsimi(0.5), val_cannot_pred(NODATA), r(0.5),
          Map_Prediction(nullptr), Map_Uncertainty(nullptr) {
    }

    ipsmIDWOperator(EnvDataset* envDataset, vector<EnvUnit *>& sampleEnvUnits)
        : CommonOperator(envDataset, sampleEnvUnits), 
		thred_envsimi(0.5), val_cannot_pred(NODATA), r(0.5),
          Map_Prediction(nullptr), Map_Uncertainty(nullptr) {
    }

    ~ipsmIDWOperator() {
    };

public:
    bool PredictMap_Property();
#ifdef MCW
    void WriteOutputs(EnvLayer* lyr, string filename);
#endif // MCW
    
public:
    double thred_envsimi;
    double val_cannot_pred;
	float r;
    EnvLayer* Map_Prediction;
    EnvLayer* Map_Uncertainty;
	string map_predict_filename;
	string map_uncer_filename;
};
}
#endif
