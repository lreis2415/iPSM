/*!
 * @brief Common operators.
 * @todo
 * @version 1.0
 * @revision  17-11-21 zhanglei - initial version
 *             17-11-27 lj       - code revision and format
 */
#ifndef COMMONOPERATOR_HPP_
#define COMMONOPERATOR_HPP_

#include <vector>
#include <math.h>
#include "DataTypeEnum.h"
#include "EnvUnit.h"
#include "EnvDataset.h"

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

namespace solim {
class CommonOperator {
public:
    EnvDataset* EDS;
    vector<EnvUnit *> SampleEnvUnits;
    vector<EnvUnit *> ValidateSampleEnvUnits;

public:
    CommonOperator() : EDS(nullptr) {
    }

    explicit CommonOperator(EnvDataset* envDataset): EDS(envDataset) {
    }

    CommonOperator(EnvDataset* envDataset, vector<EnvUnit *>& sampleEnvUnits)
     : EDS(envDataset), SampleEnvUnits(sampleEnvUnits) {
    }

    virtual ~CommonOperator();

    /*!
     * @brief 
     * @param e1
     * @param e2
     * @param mean
	 * @param stdDev
     * @param dataType
     * @return
     */
	virtual double CalcSimi_Single_Gaussian(double se, double e, double mean, double stdDev, DataTypeEnum dataType);

	/*!
	* @brief
	* @param e1
	* @param e2
	* @param range
	* @param dataType
	* @return
	*/
	virtual double CalcSimi_Single_Euclidean(double e1, double e2, double range, DataTypeEnum dataType);

    /*!
     * @brief 
     * @param e1
     * @param e2
     * @return
     */
    virtual double CalcSimi(EnvUnit* e1, EnvUnit* e2, SimilarityTypeEnum simitype = GaussianDistance);

	/*!
	* @brief
	* @param e1
	* @param e2
	* @return
	*/
	virtual double CalcSimi(EnvUnit* e1, float *envVals, SimilarityTypeEnum simitype = GaussianDistance);


    /*!
     * @brief 
     * @param e1
     * @param e2
     * @return
     */

    virtual double CalcTargetVDist(EnvUnit* e1, EnvUnit* e2);

    /*!
     * @brief 
     * @param e
     * @return
     */
    virtual double CalcUncertainty(EnvUnit* e);

	/*!
	* @brief
	* @param envVals
	* @param simitype
	* @return
	*/
	virtual double CalcUncertainty(float *envVals, SimilarityTypeEnum simitype = GaussianDistance);
};
}

#endif
