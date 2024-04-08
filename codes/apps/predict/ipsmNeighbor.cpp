#include "mpi.h"
#include "preprocess.h"
#include "ipsmNeighborOperator.h"
#include "solimIO.h"
#ifdef RasterLayer_H
#include "application.h"
using GPRO::Application;
#endif
// using namespace std; // Avoid this usage, instead of specific functions. 2019/08/06 ZHULJ
using namespace solim;

void ErrExit() {
    cout << "input format is incorrect.\n";
    exit(0);
}

int main(int argc, char *argv[]) {
    // ipsmNeighbor 参数个数：9-11
    // 输入环境因子图层文件名（以#分割），参数符号 -inlayers
    // 每个环境因子图层的数据类型，包含类别型(categorical)和连续型(continuous)，（以#分割），参数符号 -datatypes
	// 每个环境因子图层计算基于点的相似性时对应的相似性算法，包含gower、gaussian和Boolean，（以#分割），参数符号 -simimethod
    //（现在没处理这一项，默认采用高斯相似性）
    // 样点数据文件名(csv格式)，文件中包含样点地理坐标(名为:x,y,大小写均可,并且与环境因子图层的空间参考对应)和要推测的属性值，参数符号 -sample
    // 要推测的属性名(存在于样点文件的表头中)，参数符号 -target
	// 衰减系数，默认值0，参数符号 -attencoeff
    // 环境相似度阈值，默认值0.5（可选），参数符号 -simithred
    // 输出的推测图层文件名，参数符号 -predmap
    // 输出的不确定性图层文件名，参数符号 -uncmap
    // 输出各变量的特征邻域文件，参数符号 -scalemaps（可选）
    // params = "-inlayers D:/data/xc/geo.asc#D:/data/xc/planc.asc#D:/data/xc/preci.asc#D:/data/xc/profc.asc#D:/data/xc/slope.asc#D:/data/xc/tempr.asc#D:/data/xc/twi.asc -datatypes categorical#continuous#continuous#continuous#continuous#continuous#continuous -sample D:/data/xc/samples_xc.csv -target SOMB -attencoeff 2.5 -simithred 0.5 -predmap D:/data/test/pred.tif -uncmap D:/data/test/unc.tif --diffmaps";

    vector<string> envLayer_fns;
    vector<string> datatypeList;
	vector<string> similarMethod;
    string sample_fn = "./sample.csv";
	string valid_sample_fn = "./valid_sample.csv";
    string targetVName = "_";
	double attenuation_coeff = 0;
    double threshold_envSimi = 0.5;
    char* map_pred_fn = "./pred_map.tif";
    char* map_unc_fn = "./unc_map.tif";
    vector<string> scaleLayer_fns;

    //////////////////////////////////////////////////////////////////////////
    // parse params
    if (argc < 10) {
        ErrExit();
    }
    int i = 1;
    while (i < argc) {
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
        } if (strcmp(argv[i], "-simimethod") == 0) {
            i++;
            if (i < argc) {
                ParseStr(argv[i], '#', similarMethod);
                i++;
            } else { ErrExit(); }
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
        } else if (strcmp(argv[i], "-attencoeff") == 0){
			i++;
            if (i < argc) {
                sscanf(argv[i], "%lf", &attenuation_coeff);
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
        }else if (strcmp(argv[i], "-scalemaps") == 0) {
            i++;
            if (i < argc) {
                ParseStr(argv[i], '#', scaleLayer_fns);
                i++;
            } else { ErrExit(); }
        }
		else if (strcmp(argv[i], "-validsample") == 0) {
			i++;
			if (i < argc) {
				valid_sample_fn = argv[i];
				i++;
			}
			else { ErrExit(); }
		}
    }
    //////////////////////////////////////////////////////////////////////////

	// Start and initialize MPI environment
#ifdef RasterLayer_H
	Application::START(GPRO::MPI_Type, argc, argv);
	{
		int rank, size;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);


		// Begin time of loading data
		if (rank == 0) { cout << "\nLoading data ..." << endl; }
		double begint = MPI_Wtime();
#endif
		// 初始化环境因子数据
		cout << "\nloading data ...\n";
		EnvDataset *eds = new EnvDataset(envLayer_fns, datatypeList);
		EnvDataset *escales = new EnvDataset(scaleLayer_fns,datatypeList);
		vector<EnvUnit *> samples;
		vector<EnvUnit *> validPoints;
		try {
			samples = ReadTable(sample_fn, eds, targetVName);
			validPoints = ReadTable(valid_sample_fn, eds, targetVName);
		}
		catch (const char* msg) {
			std::cerr << "Input parameter error. " << msg << endl;
#ifdef RasterLayer_H
			MPI_Abort(MPI_COMM_WORLD, -1);
#else
			exit(1);
#endif
		}

#ifdef RasterLayer_H
		double readt = MPI_Wtime() - begint; // Record reading data time

											 // Create iPSM operator and perform calculation
		if (rank == 0) { cout << "\nRunning prediction method ..." << endl; }
		begint = MPI_Wtime();
#endif
		ipsmNeighborOperator *ipsmNbrOp = new ipsmNeighborOperator(eds, samples);
		ipsmNbrOp->thred_envsimi = threshold_envSimi;
		ipsmNbrOp->map_predict_filename = map_pred_fn;
		ipsmNbrOp->map_uncer_filename = map_unc_fn;


		// 调用ipsmNeighbor推理制图方法
		cout << "calculating characterisitic size layers ...";
		//ipsmNbrOp->CalcChrScale(scaleLayer_fns);
		ipsmNbrOp->Map_scale = escales;

		cout << "running prediction method ...";
		//ipsmNbrOp->PredictMap_Property(attenuation_coeff);
		ipsmNbrOp->validate(attenuation_coeff,validPoints);
		cout << "\npredicting finished!" << endl;

#ifdef RasterLayer_H
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
    cout << "\n\n--- DONE ---\n";
    return 0;
}
