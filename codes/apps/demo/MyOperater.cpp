#include "MyOperator.h"

void MyOperator::envLayer(RasterLayer<double> &layer)
{
  _pEnvLayer = &layer;
  _pEnvNbrhood = layer.nbrhood();
  cellSize = _pEnvLayer->_pMetaData->cellSize;
  noData = _pEnvLayer->_pMetaData->noData;
  Configure(_pEnvLayer, false);
}

void MyOperator::resLayer(RasterLayer<double> &layer)
{
	_pResLayer = &layer;
	Configure(_pResLayer, false);
}

bool MyOperator::Operator(const CellCoord &coord, bool operFlag)
{
	CellSpace<double> &envSpace = *(_pEnvLayer->cellSpace());
	CellSpace<double> &resSpace = *(_pResLayer->cellSpace());
	Neighborhood<double>& nbrhoodD = *(_pEnvNbrhood);
	//int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	//int dCellSize = _pEnvLayer->_pMetaData->cellSize;
	//int nodata = _pEnvLayer->_pMetaData->noData;

	int iRow = coord.iRow();
	int iCol = coord.iCol();

	int val = envSpace[iRow][iCol];

	if ( fabs(val - noData) < Eps )
	{
		resSpace[iRow][iCol] = noData;
		return true;
	}

	int res = envSpace[iRow][iCol] + 100;
	resSpace[iRow][iCol] = res;

	return true;
}
