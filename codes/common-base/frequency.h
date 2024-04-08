/*!
 * @brief
 * @version 1.0
 * @author Zhao, Fanghe
 * @revision  17-11-21 zhanglei - initial version
 */
#ifndef FREQUENCY_HPP_
#define FREQUENCY_HPP_
#include <string>
#include <list>

#include "EnvLayer.h"

// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using std::list;

namespace solim {
class Frequency {
public:
    double Value; // property value
    int Num;      // count of property value stored
    static double VERYSMALL;
public:
    Frequency() {
    }

    Frequency(double value, int num): Value(value), Num(num) {
    }

    list<Frequency>* CalFrequency(EnvLayer* lyr);

    static void PrintFrequency(list<Frequency>* freq, string filename);
};
} // namespace solim

#endif  // FREQUENCY_HPP_
