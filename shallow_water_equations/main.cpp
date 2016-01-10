#include <iostream>
#include <math.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

string to_string(int val) 
{
	char buff[32];
	sprintf(buff,"%d",val);
	return string(buff);
}

typedef struct {
	int ntimesteps; // Number of steps by time
	double dx; // Step by x
	double dy;
	double dt; // Step by time
	double g;
	int xsize; // Velocity by x
	int ysize;
	int row;
	int col;
} CommonData;

int main(int argc, char** argv)
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

	int len[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Aint pos[10] = {offsetof(CommonData, ntimesteps), offsetof(CommonData, dx), offsetof(CommonData, dy), 
					   offsetof(CommonData, dt), offsetof(CommonData, g), offsetof(CommonData, xsize),
					   offsetof(CommonData, ysize), offsetof(CommonData, row), offsetof(CommonData, col), sizeof(CommonData)};
	MPI_Datatype type[10] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_UB};

	MPI_Type_struct(10, len, pos, type, &MPI_COMMON_DATA);
	MPI_Type_commit(&MPI_COMMON_DATA);
	/* End create new MPI type */
	
	CommonData cd;
	
	double ** h;
	double startwtime, endwtime;

	if (myid == 0)
	{
		/* variables */
		cd.ntimesteps = 2000;
		cd.dx = 1;
		cd.dy = 1;
		cd.dt = 0.02;
		int xmin = 0;
		int xmax = 120;
		int ymin = 0;
		int ymax = 120;
		double x0 = 60;
		double y0 = 60;
		double h_wave = 1;
		double w_wave = -0.01;
		cd.g = 9.8;
		cd.row = 3;
		cd.col = 2;
		cd.xsize = abs(xmax - xmin)/cd.dx;
		cd.ysize = abs(ymax - ymin)/cd.dy;
	
		h = new double * [cd.xsize + 2];
		
		for (int i (0); i < cd.xsize + 2; ++i)
		{
			h[i] = new double [cd.ysize + 2];
			for (int j (0); j < cd.ysize + 2; ++j)
			{
				h[i][j] = h_wave*exp(w_wave*(pow((xmin+cd.dx*i-x0), 2) + pow((ymin+cd.dy*j-y0), 2))) + 1;
			}
		}
		system("mkdir res");
		startwtime = MPI_Wtime();
	}
	MPI_Bcast(&cd, 1, MPI_COMMON_DATA, 0, MPI_COMM_WORLD);
	if (myid != 0) 
	{
		h = new double * [cd.xsize + 2];
		for (int i (0); i < cd.xsize + 2; ++i)
		{
			h[i] = new double [cd.ysize + 2];
		}
	}
	/* variables */
	
	MPI_Bcast(&h[0][0], (cd.xsize+2)*(cd.ysize+2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double ** u = new double * [cd.xsize + 2];
	double ** v = new double * [cd.xsize + 2];
	double ** ux = new double * [cd.xsize + 1];
	double ** vx = new double * [cd.xsize + 1];
	double ** hx = new double * [cd.xsize + 1];
	double ** uy = new double * [cd.xsize + 1];
	double ** vy = new double * [cd.xsize + 1];
	double ** hy = new double * [cd.xsize + 1];
	
	for (int i (0); i < cd.xsize + 1; ++i)
	{
		hx[i] = new double [cd.ysize + 1];
		hy[i] = new double [cd.ysize + 1];
		ux[i] = new double [cd.ysize + 1];
		uy[i] = new double [cd.ysize + 1];
		vx[i] = new double [cd.ysize + 1];
		vy[i] = new double [cd.ysize + 1];
		for (int j(0); j < cd.ysize + 1; ++j)
		{
			hx[i][j] = 0;
			hy[i][j] = 0;
			ux[i][j] = 0;
			uy[i][j] = 0;
			vx[i][j] = 0;
			vy[i][j] = 0;
		}
	}

	for (int i (0); i < cd.xsize + 2; ++i)
	{
		u[i] = new double [cd.ysize + 2];
		v[i] = new double [cd.ysize + 2];
		for (int j (0); j < cd.ysize + 2; ++j)
		{
			u[i][j] = 0;
			v[i][j] = 0;
		}
	}

	ofstream out; 
	
	for (int t(0); t < cd.ntimesteps; ++t)
	{
		if (myid == 0)
		{
			for (int i(0); i < cd.xsize + 1; ++i)
			{
				for (int j (0); j < cd.ysize; ++j)
				{
					hx[i][j] = (h[i+1][j+1]+h[i][j+1])/2 - cd.dt/(2*cd.dx)*(u[i+1][j+1] - u[i][j+1]);
					ux[i][j] = (u[i+1][j+1]+u[i][j+1])/2 - cd.dt/(2*cd.dx)*((pow(u[i+1][j+1], 2)/h[i+1][j+1] + cd.g/2*pow(h[i+1][j+1], 2)) - (pow(u[i][j+1], 2)/h[i][j+1] + cd.g/2*pow(h[i][j+1], 2)));
					vx[i][j] = (v[i+1][j+1]+v[i][j+1])/2 - cd.dt/(2*cd.dx)*((u[i+1][j+1]*v[i+1][j+1]/h[i+1][j+1]) - (u[i][j+1]*v[i][j+1]/h[i][j+1]));
				}
			}
			for (int i(0); i < cd.xsize; ++i)
			{
				for (int j (0); j < cd.ysize + 1; ++j)
				{
					hy[i][j] = (h[i+1][j+1]+h[i+1][j])/2 - cd.dt/(2*cd.dy)*(v[i+1][j+1]-v[i+1][j]);
					uy[i][j] = (u[i+1][j+1]+u[i+1][j])/2 - cd.dt/(2*cd.dy)*((v[i+1][j+1]*u[i+1][j+1]/h[i+1][j+1]) - (v[i+1][j]*u[i+1][j]/h[i+1][j]));
					vy[i][j] = (v[i+1][j+1]+v[i+1][j])/2 - cd.dt/(2*cd.dy)*((pow(v[i+1][j+1], 2)/h[i+1][j+1] + cd.g/2*pow(h[i+1][j+1], 2)) - (pow(v[i+1][j], 2)/h[i+1][j] + cd.g/2*pow(h[i+1][j], 2)));
				}
			}
			for (int i(1); i < cd.xsize + 1; ++i)
			{
				for (int j(1); j < cd.ysize + 1; ++j)
				{
					h[i][j] = h[i][j] - (cd.dt/cd.dx)*(ux[i][j-1] - ux[i-1][j-1]) - (cd.dt/cd.dy)*(vy[i-1][j]-vy[i-1][j-1]);
					u[i][j] = u[i][j] - (cd.dt/cd.dx)*((pow(ux[i][j-1], 2)/hx[i][j-1] + cd.g/2*pow(hx[i][j-1], 2)) - (pow(ux[i-1][j-1], 2)/hx[i-1][j-1] + cd.g/2*pow(hx[i-1][j-1], 2)))
									  - (cd.dt/cd.dx)*((vy[i-1][j]*uy[i-1][j]/hy[i-1][j]) - (vy[i-1][j-1]*uy[i-1][j-1]/hy[i-1][j-1]));
					v[i][j] = v[i][j] - (cd.dt/cd.dx)*((ux[i][j-1]*vx[i][j-1]/hx[i][j-1]) - (ux[i-1][j-1]*vx[i-1][j-1]/hx[i-1][j-1]))
									  - (cd.dt/cd.dx)*((pow(vy[i-1][j], 2)/hy[i-1][j] + cd.g/2*pow(hy[i-1][j], 2)) - (pow(vy[i-1][j-1], 2)/hy[i-1][j-1] + cd.g/2*pow(hy[i-1][j-1], 2)));
				}
			}
		}
		else
		{
			for (int i(myid); i < ((cd.xsize + 1)/numprocs)*(myid+1); i += numprocs)
			{
				for (int j (0); j < cd.ysize; ++j)
				{
					hx[i][j] = (h[i+1][j+1]+h[i][j+1])/2 - cd.dt/(2*cd.dx)*(u[i+1][j+1] - u[i][j+1]);
					ux[i][j] = (u[i+1][j+1]+u[i][j+1])/2 - cd.dt/(2*cd.dx)*((pow(u[i+1][j+1], 2)/h[i+1][j+1] + cd.g/2*pow(h[i+1][j+1], 2)) - (pow(u[i][j+1], 2)/h[i][j+1] + cd.g/2*pow(h[i][j+1], 2)));
					vx[i][j] = (v[i+1][j+1]+v[i][j+1])/2 - cd.dt/(2*cd.dx)*((u[i+1][j+1]*v[i+1][j+1]/h[i+1][j+1]) - (u[i][j+1]*v[i][j+1]/h[i][j+1]));
				}
			}
			for (int i(myid); i < ((cd.xsize + 1)/numprocs)*(myid+1); i += numprocs)
			{
				for (int j (0); j < cd.ysize + 1; ++j)
				{
					hy[i][j] = (h[i+1][j+1]+h[i+1][j])/2 - cd.dt/(2*cd.dy)*(v[i+1][j+1]-v[i+1][j]);
					uy[i][j] = (u[i+1][j+1]+u[i+1][j])/2 - cd.dt/(2*cd.dy)*((v[i+1][j+1]*u[i+1][j+1]/h[i+1][j+1]) - (v[i+1][j]*u[i+1][j]/h[i+1][j]));
					vy[i][j] = (v[i+1][j+1]+v[i+1][j])/2 - cd.dt/(2*cd.dy)*((pow(v[i+1][j+1], 2)/h[i+1][j+1] + cd.g/2*pow(h[i+1][j+1], 2)) - (pow(v[i+1][j], 2)/h[i+1][j] + cd.g/2*pow(h[i+1][j], 2)));
				}
			}
			for (int i(myid + 1); i < ((cd.xsize + 1)/numprocs)*(myid+1); i += numprocs)
			{
				for (int j(1); j < cd.ysize + 1; ++j)
				{
					h[i][j] = h[i][j] - (cd.dt/cd.dx)*(ux[i][j-1] - ux[i-1][j-1]) - (cd.dt/cd.dy)*(vy[i-1][j]-vy[i-1][j-1]);
					u[i][j] = u[i][j] - (cd.dt/cd.dx)*((pow(ux[i][j-1], 2)/hx[i][j-1] + cd.g/2*pow(hx[i][j-1], 2)) - (pow(ux[i-1][j-1], 2)/hx[i-1][j-1] + cd.g/2*pow(hx[i-1][j-1], 2)))
									  - (cd.dt/cd.dx)*((vy[i-1][j]*uy[i-1][j]/hy[i-1][j]) - (vy[i-1][j-1]*uy[i-1][j-1]/hy[i-1][j-1]));
					v[i][j] = v[i][j] - (cd.dt/cd.dx)*((ux[i][j-1]*vx[i][j-1]/hx[i][j-1]) - (ux[i-1][j-1]*vx[i-1][j-1]/hx[i-1][j-1]))
									  - (cd.dt/cd.dx)*((pow(vy[i-1][j], 2)/hy[i-1][j] + cd.g/2*pow(hy[i-1][j], 2)) - (pow(vy[i-1][j-1], 2)/hy[i-1][j-1] + cd.g/2*pow(hy[i-1][j-1], 2)));
				}
				if (myid != 0)
				{
					MPI_Send(&h[0][0], cd.xsize*cd.ysize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
				}
			}

		}
		
		if (myid == 0)
		{
			if (numprocs > 100)
			{
				for (int i (1); i < numprocs; ++i)
				{
					double ** tmp = new double * [cd.xsize];
					for (int j(0); j < cd.xsize; ++j)
					{
						tmp[i] = new double [cd.ysize];
					}
					MPI_Recv(tmp, cd.xsize*cd.ysize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (int j(i*(cd.xsize+1)/numprocs + ((cd.xsize+1)%numprocs?1:0)); j < ((cd.xsize + 1)/numprocs + ((cd.xsize + 1)%numprocs?1:0))*(i+1); ++j)
					{
						for (int k(1); k < cd.ysize + 1; ++k)
						{
							h[j][k] = tmp[j][k];
						}
					}
				}
			}
			string name = "res/";
			name += to_string(t);
			name += ".csv";
			out.open(name.c_str());
			if (!out.is_open())
			{
				cerr << "Error! Cannot open out file!" << endl;
				return 1;
			}
			
			for (int i (0); i < cd.xsize; ++i)
			{
				for (int j(0); j < cd.ysize; ++j)
				{
					if (j == cd.ysize -1)
						out << h[i][j];
					else
						out << h[i][j] << ";";
				}
				out << "\n";
			}	

			out.close();
		}
	}
	
	for (int i(0); i < cd.xsize; ++i)
		delete h[i];
	delete h;
	
	if (myid == 0)
	{
		endwtime = MPI_Wtime();
		cout << "Done!" << endl; 
		cout << (endwtime-startwtime) << " s" <<  endl;
	}

	MPI_Finalize();

	return 0;
}
