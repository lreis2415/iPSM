#include "mpi.h"
#include "preprocess.h"
#include "ipsmFilterOperator.h"
#include "solimIO.h"

#ifdef  RasterLayer_H

#include "application.h"
using GPRO::Application;
#endif //  MPI_DEBUG

using namespace std;
using namespace solim;

void ErrExit() {
    cout << "input format is incorrect.\n";
    exit(0);
}

int main(int argc, char *argv[]) {
    // ipsm 参数个数：8
    // 输入环境因子图层文件名（以#分割），参数符号 -inlayers
    // 每个环境因子图层的数据类型，包含类别型(categorical)和连续型(continuous)，（以#分割），参数符号 -datatypes
    // 样点数据文件名(csv格式)，文件中包含样点地理坐标(名为:x,y,大小写均可,并且与环境因子图层的空间参考对应)和要推测的属性值，参数符号 -sample
    // 要推测的属性名(存在于样点文件的表头中)，参数符号 -target
    // 环境相似度阈值，默认值0.5（可选），参数符号 -simithred
    // 输出的推测图层文件名，参数符号 -predmap
    // 输出的不确定性图层文件名，参数符号 -uncmap
    // params = "-inlayers D:/data/xc/geo.asc#D:/data/xc/planc.asc#D:/data/xc/preci.asc#D:/data/xc/profc.asc#D:/data/xc/slope.asc#D:/data/xc/tempr.asc#D:/data/xc/twi.asc -datatypes categorical#continuous#continuous#continuous#continuous#continuous#continuous -sample D:/data/xc/samples_xc.csv -target SOMB -simithred 0.5 -predmap D:/data/test/pred.tif -uncmap D:/data/test/unc.tif";

    vector<string> envLayer_fns;
    vector<string> datatypeList;
	vector<string> datafactorList;
    string sample_fn = "./sample.csv";
    string targetVName = "_";
    double threshold_envSimi = 0.5;
    char* map_pred_fn = "./pred_map.tif";
    char* map_unc_fn = "./unc_map.tif";

    //////////////////////////////////////////////////////////////////////////
    // parse params
    if (argc < 9) {
        ErrExit();
    }
    int i = 1;
    while (i < argc) {
		string tmp = argv[i];
        if (strcmp(argv[i], "-inlayers") == 0) {
            i++;
            if (i < argc) {
                ParseStr(argv[i], '#', envLayer_fns);
                i++;
            } else { ErrExit(); }
        } else if (strcmp(argv[i], "-datatypes") == 0) {
            i++;
            if (i < argc) {
                ParseStr(argv[i], '#', datatypeList);
                i++;
            } else { ErrExit(); }
		}else if (strcmp(argv[i], "-datafactors") == 0) {
			i++;
			if (i < argc) {
				ParseStr(argv[i], '#', datafactorList);
				i++;
			}
			else { ErrExit(); }
        } else if (strcmp(argv[i], "-sample") == 0) {
            i++;
            if (i < argc) {
                sample_fn = argv[i];
                i++;
            } else { ErrExit(); }
        } else if (strcmp(argv[i], "-target") == 0) {
            i++;
            if (i < argc) {
                targetVName = argv[i];
                i++;
            } else { ErrExit(); }
        } else if (strcmp(argv[i], "-simithred") == 0) {
            i++;
            if (i < argc) {
                sscanf(argv[i], "%lf", &threshold_envSimi);
                i++;
            } else { ErrExit(); }
        } else if (strcmp(argv[i], "-predmap") == 0) {
            i++;
            if (i < argc) {
                map_pred_fn = argv[i];
                i++;
            } else { ErrExit(); }
        } else if (strcmp(argv[i], "-uncmap") == 0) {
            i++;
            if (i < argc) {
                map_unc_fn = argv[i];
                i++;
            } else { ErrExit(); }
		}
		else { i++; }
    }
    //////////////////////////////////////////////////////////////////////////

	// 检查相似度阈值合法性
	if (threshold_envSimi < 0.0 || threshold_envSimi > 1.0) {
		cout << "Input parameter error. Threshold of similarity must between 0 and 1.\n";
		return -1;
	}

#ifdef RasterLayer_H
	Application::START(GPRO::MPI_Type, argc, argv);
	{
#endif
#ifdef MPI_DEBUG
		int rank, size;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);


		// Begin time of loading data
		if (rank == 0) { cout << "\nLoading data ..." << endl; }
		double begint = MPI_Wtime();
#endif // MPI_DEBUG
    // 初始化环境因子数据
    cout << "\nLoading data ...\n";
    EnvDataset *eds = new EnvDataset(envLayer_fns, datatypeList);
	vector<EnvUnit *> samples;
	try {
		samples = ReadTable(sample_fn, eds, targetVName);
	}
	catch (const char* msg) {
		cerr << "Input parameter error. " << msg << endl;
#ifdef RasterLayer_H
		MPI_Abort(MPI_COMM_WORLD, -1);
#else
		exit(1);
#endif
	}
	cout << "\nLoading data ...\n";
#ifdef MPI_DEBUG
	double readt = MPI_Wtime() - begint; // Record reading data time

										 // Create iPSM operator and perform calculation
	if (rank == 0) { cout << "\nRunning prediction method ..." << endl; }
	begint = MPI_Wtime();
#endif
    ipsmFilterOperator *ipsmFilterOp = new ipsmFilterOperator(eds, samples);
	ipsmFilterOp->setDataFactors(datafactorList);
    ipsmFilterOp->thred_envsimi = threshold_envSimi;
	ipsmFilterOp->map_predict_filename = map_pred_fn;
	ipsmFilterOp->map_uncer_filename = map_unc_fn;

    // 调用ipsm推理制图方法
    cout << "\nRunning prediction method ...\n";
    ipsmFilterOp->PredictMap_Property();

    cout << "\n\n--- DONE ---\n";
#ifdef MPI_DEBUG
	double computet = MPI_Wtime() - begint; // Record computing time

											// Ouput prediction raster layers, including prediction values and uncertainty
	if (rank == 0) { cout << "\nWrite out result ..." << endl; }
	begint = MPI_Wtime();
	double writet = MPI_Wtime() - begint; // Record writing time

										  // The max. time consuming of all ranks is regarded as the actual time consuming of parallel program
	double max_readt, max_computet, max_writet;
	MPI_Allreduce(&readt, &max_readt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&computet, &max_computet, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&writet, &max_writet, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if (rank == 0) {
		printf("\nProcesses:%d\n    Read time:%f\n    Compute time:%f\n    Write time:%f\n",
			size, max_readt, max_computet, max_writet);
		fflush(stdout);
	}
#endif
#ifdef RasterLayer_H
	}
	// End and finalize MPI environment
	Application::END();
#endif
    return 0;
}