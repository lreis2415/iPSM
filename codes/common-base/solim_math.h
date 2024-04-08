/*!
 * @brief
 * @version 1.0
 * @author An Yiming
 * @revision  18-4-13
 */
#ifndef SOLIM_MATH_HPP_
#define SOLIM_MATH_HPP_
#include <vector>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::vector;

namespace solim {
class solim_math {
public:
    double envValue;
public:
    friend class ipsmNeighborOperator;

    solim_math() {
    }

    solim_math(double value) {
        envValue = value;
    }

    //get small value of two values
    int CalcMinValue(int v1, int v2);

    //calculate the similarity of two env vectors
    double CalcVectorSimi(vector<double> e1, vector<double> e2, double NoDataValue);

    //get the relative coordinate of annulus pixels to the interest pixel
    //when the outer radius of the annulus is size and the inner radius is size-1
    //the default shape of the neighborhood is circle
    void getAnnulusCord(int size, vector<int>& rCord, vector<int>& cCord);

    double CalcGaussSimi(double envPoint, double envDiff, double EnvVar);

    void CalcAnnulusWeight(double* weight, int scale, double alpha);
};
}

#endif
