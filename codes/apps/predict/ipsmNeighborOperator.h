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
#include <iostream>
#include <iomanip>
#include <numeric>

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
        void getAnnulusCord(int size, vector<int>& rCord, vector<int>& cCord);
        //calculate the similarity of two env vectors
        double CalcVectorSimi(vector<float> e1, vector<float> e2, double NoDataValue,int f);
        void getModeUnique(vector<float> v, int& mode, int& unique_num);

        // fields
    public:
        double thred_envsimi;
        double val_cannot_pred;
        int nbrMin;//min neighborhood size by anym
        int nbrMax;//max neighborhood size by anym
        EnvLayer *Map_Prediction;
        EnvLayer *Map_Uncertainty;
        string map_predict_filename;
        string map_uncer_filename;
        int n_step;
        int n_size;
        int n_step_amplifier;
        int n_range;
    };
}
#endif
