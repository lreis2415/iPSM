#ifndef _MYOPERATOR_H_
#define _MYOPERATOR_H_

#include "io_pargo.h"

using namespace GPRO;
#define Eps 0.0000001

class MyOperator : public RasterOperator<double>
{
public:
	MyOperator()
		:RasterOperator<double>(),
		_pEnvLayer(0), _pResLayer(0), num(0){}

	~MyOperator() {}

	void envLayer(RasterLayer<double> &layer);
	void resLayer(RasterLayer<double> &layer);
	virtual bool Operator(const CellCoord &coord, bool operFlag);

protected:
	double cellSize;
	double noData;
	int num;
	RasterLayer<double> *_pEnvLayer;
	RasterLayer<double> *_pResLayer;
	Neighborhood<double> *_pEnvNbrhood;
};

#endif
