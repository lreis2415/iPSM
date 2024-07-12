/*!
* @brief IPSM algorithm.
* @details �����㷨����ͷ�ļ��ļ���ͷ��ӱ�׼��ʽ��Ϣ�����������������¼��㣺by lj.
*          1. �㷨ԭ����
*          2. �ο�����
*          3. �㷨���߼�������
* @todo
* @version 1.0
* @revision  17-11-21 zhanglei - initial version
*            17-11-27 lj       - code revision and format
*/
#ifndef IPSMAIDWOPERATOR_H_
#define IPSMAIDWOPERATOR_H_

#include "CommonOperator.h"

//using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

namespace solim {
	class ipsmAIDWOperator : public CommonOperator {
	public:
		ipsmAIDWOperator() :
			thred_envsimi(0.5), val_cannot_pred(-1.),
			Map_Prediction(nullptr), Map_Uncertainty(nullptr) {
		}

		explicit ipsmAIDWOperator(EnvDataset* envDataset)
			: CommonOperator(envDataset),
			thred_envsimi(0.5), val_cannot_pred(-1.),
			Map_Prediction(nullptr), Map_Uncertainty(nullptr) {
		}

		ipsmAIDWOperator(EnvDataset* envDataset, vector<EnvUnit *>& sampleEnvUnits)
			: CommonOperator(envDataset, sampleEnvUnits),
			thred_envsimi(0.5), val_cannot_pred(-1.),
			Map_Prediction(nullptr), Map_Uncertainty(nullptr) {
		}

		~ipsmAIDWOperator() {
		};

	public:
		bool PredictMap_Property();
#ifdef MCW
		void WriteOutputs(EnvLayer* lyr, string filename);
#endif // MCW

	public:
		double thred_envsimi;
		double val_cannot_pred;
		EnvLayer* Map_Prediction;
		EnvLayer* Map_Uncertainty;
		string map_predict_filename;
		string map_uncer_filename;
	};
}
#endif
