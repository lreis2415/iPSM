#include "solim_math.h"

#include <iostream>
#include <cmath>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ

namespace solim {
//get small value of two values
int solim_math::CalcMinValue(int v1, int v2) {
    if (v1 > v2)
        return v2;
    return v1;
}

//calculate the similarity of two env vectors
double solim_math::CalcVectorSimi(vector<double> e1, vector<double> e2, double NoDataValue) {
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
void solim_math::getAnnulusCord(int size, vector<int>& rCord, vector<int>& cCord) {
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


double solim_math::CalcGaussSimi(double envPoint, double envDiff, double EnvVar) {
    double numerator = pow(envValue - envPoint, 2.0);
    double denominator = 2 * pow(EnvVar, 2.0) / envDiff;
    //cout << envValue << ", " << envPoint << ", " << envDiff <<", " << EnvVar << endl;
    double simi = exp(-numerator / denominator);
    return simi;
}



}
