/*!
 * @brief IPSMNEIGHBOR algorithm.
 * @details 对于算法，在头文件文件开头添加标准格式信息，包含但不限于以下几点：by lj.
 *          1. 算法原理简介
 *          2. 参考文献
 *          3. 算法作者及开发者
 * @todo 添加格式化代码注释. by lj..
 * @version 1.0
 * @revision  18-4-12 anyiming - initial version
 *            17-11-27 lj       - code revision and format(not done yet)
 */
#ifndef IPSMNEIGHBOROPERATOR_H_
#define IPSMNEIGHBOROPERATOR_H_

#include "CommonOperator.h"
#include "solim_math.h"//by anym
#include <iostream>
#include <iomanip>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

namespace solim
{
	class ipsmNeighborOperator : public CommonOperator
	{
		// constructor functions
	public:

		ipsmNeighborOperator() : CommonOperator()
		{
			this->Initialize();
		}

        // TODO, use explicit keyword for singe-parameter function. by lj.
		explicit ipsmNeighborOperator(EnvDataset *envDataset) : CommonOperator(envDataset)
		{
			this->Initialize();
		}

		ipsmNeighborOperator(EnvDataset *envDataset, vector<EnvUnit *> sampleEnvUnits) :
            CommonOperator(envDataset, sampleEnvUnits)
		{
			this->Initialize();
		}

    ~ipsmNeighborOperator(void) {};
		// methods
	public:
		void Initialize();

		//计算平均相似性
		//识别特征邻域大小
        bool CalcChrScale(vector<string> map_scale_names);
        bool PredictMap_Property(double alpha);
        void getEnvVector(int pCol, int pRow, EnvLayer *Layer, vector<int>rowCord,
            vector<int>colCord, vector<float> & envVector);
        //void CalcChrScale(EnvLayer * Layer);
        double CalNbrSimi_Variable(int indexPixel, int sampleRow, int sampleCol, 
			int f, float *pixel,
            EnvUnit *sample, double alpha);
        double CalNbrSimi_Sample(long pixelIndex, float *pixel, EnvUnit *sample,
            double alpha);
		bool validate(double alpha, vector<EnvUnit*>validPoints);
		
		//get the relative coordinate of annulus pixels to the interest pixel
		//when the outer radius of the annulus is size and the inner radius is size-1
		//the default shape of the neighborhood is circle
		static void getAnnulusCord(int size, vector<int>& rCord, vector<int>& cCord);
		static void CalcAnnulusWeight(double* weight, int scale, double alpha);
		//calculate the similarity of two env vectors
		static double CalcVectorSimi(vector<float> e1, vector<float> e2, double NoDataValue);



		// fields
	public:
        double thred_envsimi;
        double val_cannot_pred;
        int nbrMin;//min neighborhood size by anym
        int nbrMax;//max neighborhood size by anym
        EnvDataset *Map_scale; //anym 每个变量每个位置的特征邻域
		EnvLayer *Map_Prediction;
		EnvLayer *Map_Uncertainty;
		string map_predict_filename;
		string map_uncer_filename;
	};
}
#endif
