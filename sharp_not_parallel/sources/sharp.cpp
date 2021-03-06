#include <iostream>
#include <fstream>
#include "function.h"
#include <stdlib.h>
#include "parse_conf.h"
#include <map>

using namespace std;

void header()
{
	cout << "========================================================================" << endl;
	cout << "=======###=======###===###=======######=======########=====########=====" << endl;
	cout << "=====#######=====###===###=======######=======#########====#########====" << endl;
	cout << "===####===####===###===###======###==###======###====###===###====###===" << endl;
	cout << "===####==========###===###======###==###======###====###===###====###===" << endl;
	cout << "=====#####=======#########=====##########=====#########====#########====" << endl;
	cout << "=======#####=====#########=====##########=====########=====########=====" << endl;
	cout << "==========####===###===###====###======###====###=###======###==========" << endl;
	cout << "==========####===###===###====###======###====###==###=====###==========" << endl;
	cout << "=====#######=====###===###===###========###===###===###====###==========" << endl;
	cout << "=======###=======###===###===###========###===###====###===###==========" << endl;
	cout << "========================================================================" << endl;
	cout << "========================================================================" << endl;
	cout << "================= The software package was developed by ================" << endl;
	cout << "=================   Churaev Egor and Soldatova Elina    ================" << endl;
	cout << "=================            Copyright (c)              ================" << endl;
	cout << "========================================================================" << endl;
	cout << "========================================================================" << endl << endl;
}

int main (int argc, char ** argv)
{
	header();
	map<string, string> conf_map = parse_params(argc, argv);
	cout << "========================================================================" << endl;
	cout << "================          BEGIN PARAMETERS             =================" << endl;
	cout << "========================================================================" << endl;

	unsigned long int nx = atof(conf_map["nx"].c_str()); // Number of steps by x
	int ntimesteps = atof(conf_map["ntimesteps"].c_str()); // Number of steps by time
	double hx = atof(conf_map["hx"].c_str()); // Step by x
	double dt = atof(conf_map["dt"].c_str()); // Step by time
	double ux = atof(conf_map["ux"].c_str()); // Velocity by x
	cout << "\tnx\t\t\t\t" << nx << endl;
	cout << "\tntimesteps\t\t\t" << ntimesteps << endl;
	cout << "\thx\t\t\t\t" << hx << endl;
	cout << "\tdt\t\t\t\t" << dt << endl;
	cout << "\tux\t\t\t\t" << ux << endl;
	cout << "\tStart reading ro.init..." << endl;
	double* ro = parse_file_matrix((char*)conf_map["input_path"].c_str(), (char*)conf_map["ro"].c_str())[0];
	cout << "\tFinished reading ro.init..." << endl;
	cout << "========================================================================" << endl;
	cout << "================           END PARAMETERS              =================" << endl;
	cout << "========================================================================" << endl;
	double *fif, *dro;
	double startwtime, endwtime;

	ro = new double [nx];
	fif = new double [nx];
	dro = new double [nx];
	
	//int test_step = 2;
	//ntimesteps *= test_step;
	//dt /= test_step;

	ofstream out;
	out.open((conf_map["input_path"] + "out.csv").c_str());
	if (!out.is_open())
	{
		cerr << "Error! Cannot open out file!" << endl;
		return 1;
	}
	
	/*
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (i = myid; i <= n; i += numprocs)
	{
		drob = 1/Fact(i);
		drobSum += drob;
	}
	
	MPI_Reduce(&drobSum, &Result, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	*/
	
	cout << "========================================================================" << endl;
	cout << "================          BEGIN CALCULATING            =================" << endl;
	cout << "========================================================================" << endl;
	
	for (int t (0); t < ntimesteps; ++t)
	{
		cout << "\ttimestep\t\t\t" << t << endl;
		for (int i (2); i < nx -1; ++i)
		{
			if (ux > 0)
				fif[i] = fi_f(ro[i], ro[i-1], ro[i-2]);
			else
				fif[i] = fi_f(ro[i-1], ro[i], ro[i+1]);
		}

		ro[0] = ro[1]; //???
		for (int i (2); i < nx -1; ++i)
		{
			ro[i] += (ux * (fif[i] - fif[i + 1]) * dt/hx);
		}
		// Border conditions
		ro[nx - 1] = ro[nx - 2];
	}

	cout << "========================================================================" << endl;
	cout << "================           END CALCULATING             =================" << endl;
	cout << "========================================================================" << endl;
	
	for (int i (0); i < nx; ++i)
	{
		out << i*hx << ';' << ro[i] << '\n';
	}
	cout << "Done!" << endl; 

	out.close();
	return 0;	
}
