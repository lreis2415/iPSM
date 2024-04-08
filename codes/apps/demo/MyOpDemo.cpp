/***************************************************************************
* MyOpDemo.cpp
*
* Project: SOLIM-SOLUTION
* Purpose: A demo of basic io operations of raster files.
*
* Author:  Zhang Lei
* E-mail:  zlxy9892@163.com
****************************************************************************
* Copyright (c) 2017. Zhang Lei
*
****************************************************************************/

#include "io_pargo.h"
#include "MyOperator.h"

using namespace std;
using namespace GPRO;

int main(int argc, char *argv[])
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type,
				   Serial_Type};*/
	Application::START(MPI_Type, argc, argv); //init

	//...
	char* inputfilename = "D:/data/xc/geo.asc";
	char* neighborfile = "D:/develop/projects/SoLIM-Solutions/solim/3rdparty/pargo/nbrhoods/neigh.nbr";
	char* outputfilename = "D:/data/test/out.tif";
	//int threadNUM;

	//omp_set_num_threads(threadNUM);
	RasterLayer<double> layer1("layer1");	// ����ͼ��
	layer1.readNeighborhood(neighborfile);  // ��ȡ���������ļ�
	layer1.readFile(inputfilename);			// ��ȡդ������
	RasterLayer<double> resLayer("resLayer");
	resLayer.copyLayerInfo(layer1);

	//cout<<layer1.name();
	cout<<"\n\n";
	double starttime;
	double endtime;

	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	// TODO...
	MyOperator myOper;
	myOper.envLayer(layer1);
	myOper.resLayer(resLayer);
	//myOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	cout<<"run time is "<<endtime-starttime<<endl;

	//resLayer.writeFile(outputfilename);

	Application::END();
	return 0;
}
