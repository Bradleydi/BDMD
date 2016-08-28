//===========================================================
//===========================================================
//===========================================================
//
//  Molecular dynamics simulation of hardspheres
//
//===========================================================
//===========================================================
//===========================================================

// fa = 80 - 20   every 3 steps thermalize

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>

#include "box.h"
#include "sphere.h"
#include "event.h"
#include "heap.h"
#include "read_input.h"


int main(int argc, char **argv)
{
	read_input input;
	if (input.read(argc, argv)) return 1;

	double d, r;   // initial diameter and radius of spheres

	if(strcasecmp(input.readfile, "new")==0)
		input.readfile[0]=0;

	if (input.readfile[0]) // read in existing configuration
    {
		// read the header
		std::ifstream infile(input.readfile);
		if (!infile)
		{
			std::cout << "error, can't open " << input.readfile  << std::endl;
			exit(-1);
		}
		else
		{
			int dim;
			infile >> dim; infile.ignore(256, '\n');
			if (dim != DIM)  // quit if dimensions don't match
			{
				std::cout << "error, dimensions don't match" << std::endl;
				exit(-1);
			}
			infile.ignore(256, '\n');  // ignore the N 1 line
			infile >> input.N; infile.ignore(256, '\n');
			std::cout << "N = " << input.N << std::endl;
			infile >> d; infile.ignore(256, '\n');
			std::cout << "d = " << d << std::endl;
			r = d/2.;
			std::cout << "r = " << r << std::endl;
		}
    }
	else // create a new configuration
      r = pow(input.initialpf*pow( pow((1./6.)*M_PI*input.N/input.maxpf,1./3.)  , DIM)/(input.N*VOLUMESPHERE), 1.0/((double)(DIM)));
	
	
	box b(input.N, r, input.growthrate, input.maxpf, input.bidispersityratio, 
			input.bidispersityfraction, input.massratio, input.hardwallBC, input.polydisperse);
  
	std::cout << "ngrids = " << b.ngrids << std::endl;
	std::cout << "DIM = " << DIM << std::endl;

	if(input.readfile[0])
    {
		std::cout << "Reading in positions of spheres" << std::endl;
		b.RecreateSpheres(input.readfile, input.temp);
    }
	else 
    {
		std::cout << "Creating new positions of spheres" << std::endl;
		b.CreateSpheres(input.temp);
    } 
  
	std::ofstream output(input.datafile);
	output.precision(16);  

	b.printstep = 0;
  
	while ((b.collisionrate < input.maxcollisionrate) && (b.pf < input.maxpf) && (b.pressure < input.maxpressure)) 
    {
		b.printstep++;
		b.Process(input.N);
		//small events in 1 cycle to ensure precise final packing fraction 
		output << b.pf << " " << b.pressure << " " << 
		b.collisionrate << " " << b.neventstot << " " << std::endl;
		if(b.printstep%500==0)
			std::cout<<"time: "<<b.rtime<<" fraction: "<< b.pf <<std::endl;
		b.Synchronize(true);
    }
	output.close();

	b.WriteConfiguration(input.writefile);
	std::cout << "b.pf = " << b.pf << std::endl;
	std::cout << "b.pressure = " << b.pressure << std::endl;
	std::cout << "b.collisionrate = " << b.collisionrate << std::endl;
  
	
  
	int endtime = 50;
	int relax_time = 10;
	
	b.rtime =0;
	b.Growthrate(0.);
	b.printstep =0;
	b.Synchronize(true);
	
	b.GetRadiusAvarage();
	b.fa = 10.0;

	//FILE *stream,*flog;
	
	std::cout<<"lammpstrj file ok"<<std::endl;
	
	
	int init = 1;
	b.Process(input.N);
	
	std::cout <<"## Start EDBD relaxing... ## "<<std::endl;
	b.RUN = 1;
	
	while(b.printstep<5000)//relax
	{
		b.printstep++;
		//std::cout <<printstep<<std::endl;
		//b.Process(input.N);
		if(init)
		{
			b.Synchronize_assign();//including check and find event
			b.Process(input.N);
			init--;
		}
		
			if(b.printstep%2==0)
			{
				
				b.Synchronize_update_fa(); //including check and find event
				b.Process(input.N/2);//run

				if(b.printstep%60==0)
				{
					std::cout<<"time: "<<b.rtime<<" fraction: "<< b.pf <<std::endl;
				}
			}
			
			else
			{
				b.ncollisions = 0;	b.neventstot = 0;
				b.Process(input.N/2);
				//std::cout<<"collision "<<b.ncollisions<<" collision tot "<<b.neventstot<<std::endl;
			}
		
		//b.Synchronize(true);
	}
	
	std::cout<<"deltat: "<<b.deltat2<<std::endl;
	
	double dtemp=10.;
	while(floor(b.deltat2*dtemp)==0)dtemp *= 10.;
	b.deltat = floor(b.deltat2*dtemp)/dtemp;
	
	std::cout<<"We choose deltat: "<<b.deltat<<std::endl;
	std::cout<<"## Start EDBD ##"<<std::endl;
	
	
	b.RUN_2 = 1;
	b.rtime = 0;	 b.time_old = 0;	b.printstep = 0;
	b.btime = b.rtime + b.deltat;
	//b.Process(input.N);
	b.Assign_xu();
	
	init=1;
	while(b.rtime < 10000.0)
	{
		if(init){
			b.PrintLog();
			b.PrintProcess(b.printstep);
			init--;
		}
		//std::cout <<printstep<<std::endl;
		//output << b.pf << " " << b.pressure << " " << 
		//b.collisionrate << " " << b.neventstot << " " << std::endl;
		b.Process(input.N);	
		
	
	}
  fclose(b.flog);
  fclose(b.stream);
  return 0;
  
}
