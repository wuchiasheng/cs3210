//
//  metropolis.c
//  
//
//  Created by Wu Chia Sheng on 21/4/14.
//
//

//particle is initialized such that the distance is greater than 1 hence attractive energy, relative fast to find lesser energy
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define N 32

typedef struct {
	double x;
	double y;
	double z;
}particle;

void createParticles(particle p[N]){
	int i;
	
	/* initialize random seed: */
  	srand (time(NULL));
	
    
	for(i=0;i<N;i++){
		p[i].x=(rand() % 100 + 1)/10.00;
		p[i].y=(rand() % 100 + 1)/10.00;
		p[i].z=(rand() % 100 + 1)/10.00;
	}
    
}

int main(int argc, char** argv){
	
	int i,j,k,l;
	
	int timeStep=0;
	int worldRank;
	int worldSize;
	
	double elapseTime;
	clock_t start;
	clock_t end;
	float seconds1=0.00;
	float seconds3=0.00;
	
	int numParticle;
	int dispParticle;
	double distance;
	double initialEnergy;
	double energy;
	double totalEnergy=0.00;
	double localEnergy=0.00;
	double globalEnergy;
	double previousEnergy=0.00;
	double fabGlobal;
	double fabPrev;
	particle p[N];
	particle oldP[N];
	int newX,newY,newZ;
	
	//496 is derived from 32P2/2
	int index1Pair[496];
	int index2Pair[496];
	
	MPI_Datatype particleType;
	MPI_Datatype type[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	int blocklen[3]={1,1,1};
	MPI_Aint disp[3];
	MPI_Aint extent;
	//Initialize MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
	MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
	MPI_Type_extent(MPI_DOUBLE,&extent);
	disp[0]=0;
	disp[1]=1*extent;
	disp[2]=2*extent;
	MPI_Type_struct(3,blocklen,disp,type,&particleType);
	MPI_Type_commit(&particleType);
    
	if(worldRank==0){
		createParticles(p);
		for(i=0;i<N;i++)
		{
			printf("Particle %d: X:%f,Y:%f,Z:%f\n",i,p[i].x,p[i].y,p[i].z);
		}
	}
	
	while(1){
		totalEnergy=0.00;
		localEnergy=0.00;
		globalEnergy=0.00;
		//broadcast particle coordinates to all processors
		MPI_Bcast(p,N,particleType, 0, MPI_COMM_WORLD);
        
		/*assign  number of particle to compute
         particle decomposition technique is used
         */
		numParticle =31;
		dispParticle = worldRank*(numParticle);
        
		//calculate total energy
        
		//form array 1 that contains the index of first particle of a pair
		l=0;
		k=N-1;//-1 to exclude same index pair
        
		for(i=0;i<N-1;i++){
			for(j=0;j<k;j++){
				index1Pair[l]=(N-1)-k+1;
				l++;
			}
			k--;
		}
        
		
        //form array 2 that contains the index of second particle of a pair
		l=0;
		k=N-1;
		for(i=0;i<N-1;i++){
			for(j=0;j<k;j++){
				index2Pair[l]=(N-1)-k+2+j;
				l++;
			}
			k--;
		}
   
		for(i=0;i<numParticle;i++){
			k=i+dispParticle;
			distance = sqrt(pow(p[index1Pair[k]].x-p[index2Pair[k]].x,2.0)+pow(p[index1Pair[k]].y-p[index2Pair[k]].y,2.0)+pow(p[index1Pair[k]].z-p[index2Pair[k]].z,2.0));
			if(distance == 0){
				energy=0.00;
			}
			else
				energy = 1/pow(distance,12.0)-1/pow(distance,6.0);
			totalEnergy = totalEnergy + energy;
		}
		localEnergy = totalEnergy;
        
		//sum of all energy
		MPI_Allreduce(&localEnergy,&globalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        
		if(timeStep ==0){
			initialEnergy=globalEnergy;
		}
		if(timeStep ==0 && worldRank==0)
			printf("Initial Energy= %f\n",initialEnergy);
        
		if(timeStep ==0){
            initialEnergy=fabs(globalEnergy);
        }
	
		
		//accept new coordinate if new state energy less than previous state
		fabGlobal = fabs(globalEnergy);
		fabPrev = fabs(previousEnergy);
        
		if( fabGlobal < fabPrev || timeStep==50000){
			if(fabGlobal <= (initialEnergy/4.00) && timeStep!=0 || timeStep ==50000){
				if(worldRank==0){
                    printf("Solution is found after %d\n",timeStep);
                    printf("New Energy = %f\n",globalEnergy);
				}
				break;
			}
		}
		
		else if(timeStep!=0){
			//if(worldRank==0)
            //printf("greater");
            for(i=0;i<N;i++){
                p[i].x=oldP[i].x;
                p[i].y=oldP[i].y;
                p[i].z=oldP[i].z;
                
            }
		}
        
		previousEnergy = globalEnergy;
        
		//save old coordinates before moving particles
		for(i=0;i<N;i++){
			oldP[i].x=p[i].x;
			oldP[i].y=p[i].y;
			oldP[i].z=p[i].z;
            
		}
        
		/*
         MOVING PARTICLES
         */
		if(worldRank==0){
			srand (time(NULL));
			for(i=0;i<N;i++){
				int randX=rand() % 19 + (-9);
				int randY=rand() % 19 + (-9);
				int randZ=rand() % 19 + (-9);
				newX = p[i].x*100 + randX*1;
				newY = p[i].y*100 + randY*1;
				newZ = p[i].z*100 + randZ*1;
				
				if(newX<=900 && newX>=100)
					p[i].x = (double)newX/100;
				if(newY<=900 && newY>=100)
					p[i].y = (double)newY/100;
				if(newZ<=900 && newZ>=100)
					p[i].z = (double)newZ/100;
			}
            
         
		}
		timeStep++;		
        
	}
    
	MPI_Finalize();
	return 0;
}
