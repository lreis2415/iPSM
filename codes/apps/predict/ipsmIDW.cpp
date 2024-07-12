#include "mpi.h"

#include "preprocess.h"
#include "solimIO.h"
#include "ipsmIDWOperator.h"
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
    // params = "-inlayers D:/data/xc/geo.asc#D:/data/xc/planc.asc#D:/data/xc/preci.asc#D:/data/xc/profc.asc#D:/data/xc/slope.asc#D:/data/xc/tempr.asc#D:/data/xc/twi.asc -datatypes categorical#continuous#continuous#continuous#continuous#continuous#continuous -sample D:/data/xc/samples_xc.csv -target SOMB -simithred 0.5 -predmap D:/data/test/pred.tif -uncmap D:/data/test/unc.tif";

    vector<string> envLayer_fns;
    vector<string> datatypeList;
    string sample_fn = "./sample.csv";
    string targetVName = "_";
    double threshold_envSimi = 0.5;
	float rValue = 0.75;
    char* map_pred_fn = "./pred_map.tif";
    char* map_unc_fn = "./unc_map.tif";

    //////////////////////////////////////////////////////////////////////////
    // parse params
    if (argc < 8) {
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
    }
    //////////////////////////////////////////////////////////////////////////

    if (threshold_envSimi < 0.0 || threshold_envSimi > 1.0) {
        cout << "Input parameter error. Threshold of similarity must between 0 and 1.\n";
        return -1;
    }
    // Start and initialize MPI environment
#ifdef RasterLayer_H
	Application::START(GPRO::MPI_Type, argc, argv);
	{
#endif
#ifdef RasterLayer_H
		int rank, size;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);


		// Begin time of loading data
		if (rank == 0) { cout << "\nLoading data ..." << endl; }
		double begint = MPI_Wtime();
#endif // RasterLayer_H
        // Initialize environment variables from GeoTiff
        EnvDataset *eds = new EnvDataset(envLayer_fns, datatypeList);
        // Read samples from table (e.g., .csv file)
        vector<EnvUnit *> samples;
        try {
            samples = ReadTable(sample_fn, eds, targetVName);
			
        } catch (const char* msg) {
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
		//cout << "sample size:" << samples.size() << endl;
        ipsmIDWOperator *ipsmOp = new ipsmIDWOperator(eds, samples);
        ipsmOp->thred_envsimi = threshold_envSimi;
		ipsmOp->r = rValue;
		cout << "pre:" << map_pred_fn;
		cout << "unc:" << map_unc_fn << endl;
		ipsmOp->map_predict_filename = map_pred_fn;
		ipsmOp->map_uncer_filename = map_unc_fn;
        ipsmOp->PredictMap_Property();
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
	return 0;
}
