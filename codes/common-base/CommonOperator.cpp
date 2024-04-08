#include "CommonOperator.h"

#include <cassert>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

namespace solim {

CommonOperator::~CommonOperator() {
    if (EDS != nullptr) {
        delete EDS;
        EDS = nullptr;
    }
}

double CommonOperator::CalcSimi(EnvUnit* se, EnvUnit* e, SimilarityTypeEnum simitype) {
    double simi = -1.;
    if (se->EnvValues.size() != e->EnvValues.size()) return -1.;
    if (!se->IsCal || !e->IsCal) return -1.;
	simi = 1;
	if (simitype == GaussianDistance) {
		for (int i = 0; i < se->EnvValues.size(); i++) {
			//double range = EDS->Layers.at(i)->Data_Range;
			double mean = EDS->Layers.at(i)->Data_Mean;
			double stdDev = EDS->Layers.at(i)->Data_StdDev;
			DataTypeEnum dataType = se->DataTypes.at(i);
			double simi_temp = CalcSimi_Single_Gaussian(se->EnvValues.at(i), e->EnvValues.at(i), mean, stdDev, dataType);
			if (simi_temp < simi) {
				simi = simi_temp;
			}
		}
    }
	else if (simitype == EuclideanDistance) {
		for (int i = 0; i < se->EnvValues.size(); i++) {
			double range = EDS->Layers.at(i)->Data_Range;
			DataTypeEnum dataType = se->DataTypes.at(i);
			double simi_temp = CalcSimi_Single_Euclidean(se->EnvValues.at(i), e->EnvValues.at(i), range, dataType);
			if (simi_temp < simi) {
				simi = simi_temp;
			}
		}
	}
    return simi;
}

double CommonOperator::CalcSimi(EnvUnit* se, float* envVals, SimilarityTypeEnum simitype) {
	double simi = -1.;
	if (se->EnvValues.size() != EDS->LayerSize) return -1.;
	if (!se->IsCal) return -1.;
	for (int i = 0; i < se->EnvValues.size(); i++) {
		float nodata = EDS->Layers.at(i)->NoDataValue;
		if (fabs(envVals[i] - nodata) < VERY_SMALL /*|| isnan(envVals[i] - nodata)*/) return -1;
	}
	simi = 1;
	if (simitype == GaussianDistance) {
		for (int i = 0; i < se->EnvValues.size(); i++) {
			double mean = EDS->Layers.at(i)->Data_Mean;
			double stdDev = EDS->Layers.at(i)->Data_StdDev;
			DataTypeEnum dataType = se->DataTypes.at(i);
			double simi_temp = CalcSimi_Single_Gaussian(se->EnvValues.at(i), envVals[i], mean, stdDev, dataType);
			if (simi_temp < simi) {
				simi = simi_temp;
			}
		}
	}
	else if (simitype == EuclideanDistance) {
		for (int i = 0; i < se->EnvValues.size(); i++) {
			double range = EDS->Layers.at(i)->Data_Range;
			DataTypeEnum dataType = se->DataTypes.at(i);
			double simi_temp = CalcSimi_Single_Euclidean(se->EnvValues.at(i), envVals[i], range, dataType);
			if (simi_temp < simi) {
				simi = simi_temp;
			}
		}
	}
	return simi;
	
	return -1.;
}


double CommonOperator::CalcSimi_Single_Gaussian(double se, double e, double mean, double stdDev, DataTypeEnum dataType) {
    //assert(range != 0.);

    double simi = -1.;
    if (dataType == CATEGORICAL) {
        if (se == e) {
            simi = 1.;
        } else {
            simi = 0.;
        }
    } else if (dataType == CONTINUOUS) {
        //simi = 1. - abs(se - e) / range;
		simi = exp(-pow(se - e, 2) * 0.5 / pow(stdDev, 4) * (se * se + stdDev * stdDev + mean * mean - 2 * se*mean));
    }
    return simi;
}

double CommonOperator::CalcSimi_Single_Euclidean(double e1, double e2, double range, DataTypeEnum dataType) {
	assert(range != 0.);

	double simi = -1.;
	if (dataType == CATEGORICAL) {
		if (e1 == e2) {
			simi = 1.;
		}
		else {
			simi = 0.;
		}
	}
	else if (dataType == CONTINUOUS) {
		simi = 1. - abs(e1 - e2) / range; //membership function is linear relationship
	}
	return simi;
}

double CommonOperator::CalcTargetVDist(EnvUnit* e1, EnvUnit* e2) {
    double dist = 0;
    dist = abs(e1->SoilVariable - e2->SoilVariable);
    return dist;
}

double CommonOperator::CalcUncertainty(EnvUnit* e) {
    double simi = 0.;
    double simi_temp = 0.;
    for (vector<EnvUnit *>::iterator it = SampleEnvUnits.begin(); it != SampleEnvUnits.end(); ++it) {
        simi_temp = CalcSimi(*it, e);
        if (simi < simi_temp) {
            simi = simi_temp;
        }
    }
    return 1. - simi;
}

double CommonOperator::CalcUncertainty(float* envVals, SimilarityTypeEnum simitype) {
	double simi = 0.;
	double simi_temp = 0.;
	for (vector<EnvUnit *>::iterator it = SampleEnvUnits.begin(); it != SampleEnvUnits.end(); ++it) {
		simi_temp = CalcSimi(*it, envVals, simitype);
		if (simi < simi_temp) {
			simi = simi_temp;
		}
	}
	return 1. - simi;
}
}

