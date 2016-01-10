#include <iostream>
#include <fstream>
#include "function.h"
#include <stdlib.h>
#include "parse_conf.h"
#include <map>
#include <mpi.h>

using namespace std;

typedef struct {
	unsigned long int nx; // Number of steps by x
	int ntimesteps; // Number of steps by time
	double hx; // Step by x
	double dt; // Step by time
	double ux; // Velocity by x
} CommonData;

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
	int init;
	int myid;
	int numprocs;
	if (init = MPI_Init(&argc, &argv))
	{
		cerr << "Error! Cannot run MPI!" << endl;
		MPI_Abort(MPI_COMM_WORLD, init);
	}
	
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* Start create new MPI type */
	MPI_Datatype MPI_COMMON_DATA;

	int len[6] = {1, 1, 1, 1, 1, 1};
	MPI_Aint pos[6] = {offsetof(CommonData, nx), offsetof(CommonData, ntimesteps), offsetof(CommonData, hx), 
					   offsetof(CommonData, dt), offsetof(CommonData, ux), sizeof(CommonData)};
	MPI_Datatype type[6] = {MPI_UNSIGNED_LONG, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UB};

	MPI_Type_struct(6, len, pos, type, &MPI_COMMON_DATA);
	MPI_Type_commit(&MPI_COMMON_DATA);
	/* End create new MPI type */
	
	CommonData cd;
	double* ro;
	double startwtime, endwtime;
	ofstream out;
	int ro_num = 0;

	if (myid == 0) 
	{
		header();
		map<string, string> conf_map = parse_params(argc, argv);
		cout << "========================================================================" << endl;
		cout << "================          BEGIN PARAMETERS             =================" << endl;
		cout << "========================================================================" << endl;

		cd.nx = atof(conf_map["nx"].c_str()); // Number of steps by x
		cd.ntimesteps = atof(conf_map["ntimesteps"].c_str()); // Number of steps by time
		cd.hx = atof(conf_map["hx"].c_str()); // Step by x
		cd.dt = atof(conf_map["dt"].c_str()); // Step by time
		cd.ux = atof(conf_map["ux"].c_str()); // Velocity by x
		cout << "\tnx\t\t\t\t" << cd.nx << endl;
		cout << "\tntimesteps\t\t\t" << cd.ntimesteps << endl;
		cout << "\thx\t\t\t\t" << cd.hx << endl;
		cout << "\tdt\t\t\t\t" << cd.dt << endl;
		cout << "\tux\t\t\t\t" << cd.ux << endl;
		cout << "\tStart reading ro.init..." << endl;
		ro = parse_file_matrix((char*)conf_map["input_path"].c_str(), (char*)conf_map["ro"].c_str())[0];
		ro_num = cd.nx;
		cout << "\tFinished reading ro.init..." << endl;
		cout << "========================================================================" << endl;
		cout << "================           END PARAMETERS              =================" << endl;
		cout << "========================================================================" << endl;
		out.open((conf_map["input_path"] + "out.csv").c_str());
		if (!out.is_open())
		{
			cerr << "Error! Cannot open out file!" << endl;
			return 1;
		}

		cout << "========================================================================" << endl;
		cout << "================          BEGIN CALCULATING            =================" << endl;
		cout << "========================================================================" << endl;
		startwtime = MPI_Wtime();
	}
	
	// Send init info to all process
	MPI_Bcast(&cd, 1, MPI_COMMON_DATA, 0, MPI_COMM_WORLD);
	
	if (myid != 0)
	{
		ro = new double [cd.nx];
	}
	double* fif = new double [cd.nx];
	const int buf_size = cd.nx / numprocs + ((cd.nx % numprocs)?1:0);
	double* recv_buf = new double [buf_size];
	double* send_buf = new double [buf_size];
	int send_to, recv_from;
	if (myid + 1 == numprocs)
		send_to = 0;
	else
		send_to = myid + 1;

	if (myid - 1 < 0)
		recv_from = numprocs - 1;
	else
		recv_from = myid - 1;

	bool final_proc = false;
	for (int t(myid); t < cd.ntimesteps; t += numprocs)
	{
		int sent_ro = 0;
		int start_ro = 2;
		if (numprocs == 1)
		{
			ro_num = cd.nx;
		}
		do
		{
			if (ro_num + buf_size < cd.nx)
			{
				MPI_Recv(recv_buf, buf_size, MPI_DOUBLE, recv_from, t-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i(0); i < buf_size; ++i)
				{
					ro[ro_num + i] = recv_buf[i];
				}
				ro_num += buf_size;
			}
			else if (ro_num < cd.nx)
			{
				MPI_Recv(recv_buf, cd.nx-ro_num, MPI_DOUBLE, recv_from, t-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i(0); i < cd.nx - ro_num; ++i)
				{
					ro[ro_num + i] = recv_buf[i];
				}
				ro_num += (cd.nx - ro_num);
			}
			for (int i (start_ro); i < ro_num-1; ++i)
			{
				if (cd.ux > 0)
					fif[i] = fi_f(ro[i], ro[i-1], ro[i-2]);
				else
					fif[i] = fi_f(ro[i-1], ro[i], ro[i+1]);
			}

			if (start_ro == 2)
			{
				ro[0] = ro[1]; 
			}
			for (int i (start_ro); i < ro_num -2; ++i)
			{
				ro[i] += (cd.ux * (fif[i] - fif[i + 1]) * cd.dt/cd.hx);
				if (numprocs > 1 && t < cd.ntimesteps - 1) 
				{
					if (i % buf_size == 0)
					{
						for (int j(0); j < buf_size; ++j)
						{
							send_buf[j] = ro[i - buf_size + j];
						}
						MPI_Send(send_buf, buf_size, MPI_DOUBLE, send_to, t, MPI_COMM_WORLD);
						sent_ro++;
					}
				}
			}
			if (ro_num == cd.nx)
			{
				ro[cd.nx - 2] += (cd.ux * (fif[cd.nx - 2] - fif[cd.nx - 1]) * cd.dt/cd.hx);
				// Border conditions
				ro[cd.nx - 1] = ro[cd.nx - 2];
				if (numprocs > 1 && t < cd.ntimesteps - 1)
				{
					for (int j(0); j < cd.nx - sent_ro*buf_size; ++j)
					{
						send_buf[j] = ro[sent_ro*buf_size + j];
					}
					MPI_Send(send_buf, cd.nx - sent_ro*buf_size, MPI_DOUBLE, send_to, t, MPI_COMM_WORLD);
				}
			}
			start_ro = ro_num;
		}
		while (ro_num < cd.nx);
		cout << "\ttimestep\t\t\t" << t << endl;
		ro_num = 0;
		if (t == cd.ntimesteps - 1)
		{
			final_proc = true;
		}
	}
	
	delete recv_buf;
	delete send_buf;
	delete fif;

	if (final_proc && myid != 0)
	{
		MPI_Send(ro, cd.nx, MPI_DOUBLE, 0, 888, MPI_COMM_WORLD);
	}

	if (myid == 0)
	{   
		if (!final_proc)
		{
			MPI_Recv(ro, cd.nx, MPI_DOUBLE, MPI_ANY_SOURCE, 888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		cout << "========================================================================" << endl;
		cout << "================           END CALCULATING             =================" << endl;
		cout << "========================================================================" << endl;
		endwtime = MPI_Wtime();
		for (int i (0); i < cd.nx; ++i)
		{
			out << i*cd.hx << ';' << ro[i] << '\n';
		}
		out.close();
		cout << "Done!" << endl; 
		cout << (endwtime-startwtime) << " s" <<  endl;
	}
	delete ro;
	MPI_Finalize();
}
