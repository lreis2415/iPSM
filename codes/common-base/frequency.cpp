#include "frequency.h"

#include <iostream>
#include <string>
#include <list>
#include <fstream>
#include <cmath>

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::list;
using std::ofstream;

namespace solim {
double Frequency::VERYSMALL = 0.0001;
list<Frequency>* Frequency::CalFrequency(EnvLayer* lyr) {
    // 计算样本属性值出现的频次
    list<Frequency>* freq = new list<Frequency>;
    double val;
    int first_loc = 0;
    while (1) {
        if (fabs(lyr->EnvData[first_loc] - lyr->NoDataValue) < VERYSMALL) {
            first_loc++;
        } else {
            break;
        }
    }
    freq->push_back(Frequency(lyr->EnvData[first_loc], 1));
    for (int i = first_loc + 1; i < lyr->XSize * lyr->YSize; i++) {
        val = lyr->EnvData[i];
        if (abs(val - lyr->NoDataValue) < VERYSMALL) { continue; }
        bool greatest = true;
        for (list<Frequency>::iterator itr = freq->begin(); itr != freq->end(); ++itr) {
            if (abs(itr->Value - val) < VERYSMALL) {
                //If the Value of current iterator equals to val(current EnvData value), add 1 to Num of the current iterator
                itr->Num++;
                greatest = false;
                break;
            }

            if (val < itr->Value) {
                //If the val is smaller than the Value of current iterator, insert Frequency(val,1) before current iterator
                freq->insert(itr, Frequency(val, 1));
                greatest = false;
                break;
            }
            //If the val is bigger than the Value of current iterator, circulation goes on
        }
        if (greatest) {
            //如果当前属性值最大，将其插入表尾
            freq->push_back(Frequency(val, 1));
        }
    }
    return freq;
}

void Frequency::PrintFrequency(list<Frequency>* freq, string filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cout << "Cannot open output file" << endl;
        return;
    }
    for (list<Frequency>::iterator itr = freq->begin(); itr != freq->end(); itr++) {
        outfile << itr->Value << " " << itr->Num << endl;
    }
    freq->clear();
    freq = NULL;
}
};
