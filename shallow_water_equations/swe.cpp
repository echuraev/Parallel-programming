#include <iostream>
#include <math.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>

using namespace std;

string to_string(int val) 
{
	char buff[32];
	sprintf(buff,"%d",val);
	return string(buff);
}

int main()
{
	/* variables */
	int ntimestep = 2000;
	double dx = 1;
	double dy = 1;
	double dt = 0.02;
	int xmin = 0;
	int xmax = 120;
	int ymin = 0;
	int ymax = 120;
	double x0 = 60;
	double y0 = 60;
	double h_wave = 1;
	double w_wave = -0.01;
	double g = 9.8;
	int xsize = abs(xmax - xmin)/dx;
	int ysize = abs(ymax - ymin)/dy;
	
	double ** h = new double * [xsize + 2];
	double ** u = new double * [xsize + 2];
	double ** v = new double * [xsize + 2];
	double ** ux = new double * [xsize + 1];
	double ** vx = new double * [xsize + 1];
	double ** hx = new double * [xsize + 1];
	double ** uy = new double * [xsize + 1];
	double ** vy = new double * [xsize + 1];
	double ** hy = new double * [xsize + 1];
	
	for (int i (0); i < xsize + 1; ++i)
	{
		hx[i] = new double [ysize + 1];
		hy[i] = new double [ysize + 1];
		ux[i] = new double [ysize + 1];
		uy[i] = new double [ysize + 1];
		vx[i] = new double [ysize + 1];
		vy[i] = new double [ysize + 1];
		for (int j(0); j < ysize + 1; ++j)
		{
			hx[i][j] = 0;
			hy[i][j] = 0;
			ux[i][j] = 0;
			uy[i][j] = 0;
			vx[i][j] = 0;
			vy[i][j] = 0;
		}
	}
	for (int i (0); i < xsize + 2; ++i)
	{
		h[i] = new double [ysize + 2];
		u[i] = new double [ysize + 2];
		v[i] = new double [ysize + 2];
		for (int j (0); j < ysize + 2; ++j)
		{
			h[i][j] = h_wave*exp(w_wave*(pow((xmin+dx*i-x0), 2) + pow((ymin+dy*j-y0), 2))) + 1;
			u[i][j] = 0;
			v[i][j] = 0;
		}
	}
	/* variables */
	
	ofstream out; 
	
	system("mkdir res");
	for (int t(0); t < ntimestep; ++t)
	{
		for (int i(0); i < xsize + 1; ++i)
		{
			for (int j (0); j < ysize; ++j)
			{
				hx[i][j] = (h[i+1][j+1]+h[i][j+1])/2 - dt/(2*dx)*(u[i+1][j+1] - u[i][j+1]);
				ux[i][j] = (u[i+1][j+1]+u[i][j+1])/2 - dt/(2*dx)*((pow(u[i+1][j+1], 2)/h[i+1][j+1] + g/2*pow(h[i+1][j+1], 2)) - (pow(u[i][j+1], 2)/h[i][j+1] + g/2*pow(h[i][j+1], 2)));
				vx[i][j] = (v[i+1][j+1]+v[i][j+1])/2 - dt/(2*dx)*((u[i+1][j+1]*v[i+1][j+1]/h[i+1][j+1]) - (u[i][j+1]*v[i][j+1]/h[i][j+1]));
			}
		}
		for (int i(0); i < xsize; ++i)
		{
			for (int j (0); j < ysize + 1; ++j)
			{
				hy[i][j] = (h[i+1][j+1]+h[i+1][j])/2 - dt/(2*dy)*(v[i+1][j+1]-v[i+1][j]);
				uy[i][j] = (u[i+1][j+1]+u[i+1][j])/2 - dt/(2*dy)*((v[i+1][j+1]*u[i+1][j+1]/h[i+1][j+1]) - (v[i+1][j]*u[i+1][j]/h[i+1][j]));
				vy[i][j] = (v[i+1][j+1]+v[i+1][j])/2 - dt/(2*dy)*((pow(v[i+1][j+1], 2)/h[i+1][j+1] + g/2*pow(h[i+1][j+1], 2)) - (pow(v[i+1][j], 2)/h[i+1][j] + g/2*pow(h[i+1][j], 2)));
			}
		}
		for (int i(1); i < xsize + 1; ++i)
		{
			for (int j(1); j < ysize + 1; ++j)
			{
				//h[i][j] = (1 - omega) * h[i][j] + omega * (h[i][j+1] + h[i][j-1] + h[i+1][j] + h[i-1][j])/4.0;
				h[i][j] = h[i][j] - (dt/dx)*(ux[i][j-1] - ux[i-1][j-1]) - (dt/dy)*(vy[i-1][j]-vy[i-1][j-1]);
				u[i][j] = u[i][j] - (dt/dx)*((pow(ux[i][j-1], 2)/hx[i][j-1] + g/2*pow(hx[i][j-1], 2)) - (pow(ux[i-1][j-1], 2)/hx[i-1][j-1] + g/2*pow(hx[i-1][j-1], 2)))
								  - (dt/dx)*((vy[i-1][j]*uy[i-1][j]/hy[i-1][j]) - (vy[i-1][j-1]*uy[i-1][j-1]/hy[i-1][j-1]));
				v[i][j] = v[i][j] - (dt/dx)*((ux[i][j-1]*vx[i][j-1]/hx[i][j-1]) - (ux[i-1][j-1]*vx[i-1][j-1]/hx[i-1][j-1]))
								  - (dt/dx)*((pow(vy[i-1][j], 2)/hy[i-1][j] + g/2*pow(hy[i-1][j], 2)) - (pow(vy[i-1][j-1], 2)/hy[i-1][j-1] + g/2*pow(hy[i-1][j-1], 2)));
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
		
		for (int i (0); i < xsize; ++i)
		{
			for (int j(0); j < ysize; ++j)
			{
				if (j == ysize -1)
					out << h[i][j];
				else
					out << h[i][j] << ";";
			}
			out << "\n";
		}

		out.close();
	}
	
	for (int i(0); i < xsize; ++i)
		delete h[i];
	delete h;
	
	return 0;
}
