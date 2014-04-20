//
//  energycalculation.c
//  
//
//  Created by Wu Chia Sheng on 21/4/14.
//
//

//this piece of codes implement 2 algorithms to calculate total energy of a system
//this piece of codes also reflect the timing of the 2 algorithms for comparison

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
		p[i].x=(rand() % 10 + 1);
		p[i].y=(rand() % 10 + 1);
		p[i].z=(rand() % 10 + 1);
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
	particle p[N];
	particle oldP[N];
	int newX,newY,newZ;
	
	//496 is derived by 32P2/2
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
    
	int step=0;
    
	if(worldRank==0){
		createParticles(p);
		for(i=0;i<N;i++)
		{
			printf("Particle %d: X:%f,Y:%f,Z:%f\n",i,p[i].x,p[i].y,p[i].z);
		}
	}
	
	while(step<1000){
        step++;
        MPI_Bcast(p,N,particleType, 0, MPI_COMM_WORLD);
        
        //-----------------------------------------------------algorithm 2 energy calculation-----------------------------
        /*assign  number of particle to compute
         particle decomposition technique is used
         */
        numParticle =31;
        dispParticle = worldRank*(numParticle);
        
        //calculate total energy
        start = clock();
        
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
        
        //starts calculating energy between 2 particles assigned
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
        end = clock();
        seconds3 =seconds3 + (float) (end-start) / CLOCKS_PER_SEC;
        
        //sum of all energy
        MPI_Allreduce(&localEnergy,&globalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if(worldRank==0 && step==1)
            printf("Total energy from algo3 %d = %f\n",worldRank,globalEnergy);
        
        //----------------------------------------Algorithm 1 energy calculation-------------------------------------------
        
        //calculate number of particle to compute
        numParticle =2;
        dispParticle = worldRank*(numParticle);
        localEnergy=0.00;
        globalEnergy=0.00;
        
        
        start = clock();
        
        for(i=0;i<numParticle;i++){
            totalEnergy=0.00;
            k=i+dispParticle;
            for(j=0;j<N;j++){
                if(j!=k){
                    distance = sqrt(pow(p[k].x-p[j].x,2.0)+pow(p[k].y-p[j].y,2.0)+pow(p[k].z-p[j].z,2.0));
                    if(distance == 0){
                        energy=0.00;
                    }
                    else
                        energy = 1/pow(distance,12.0)-1/pow(distance,6.0);
                    totalEnergy = totalEnergy + energy;
                }
                
            }
            localEnergy =localEnergy + totalEnergy;
        }
        
        end = clock();
        seconds1 = seconds1 + (float) (end-start) / CLOCKS_PER_SEC;
        
        
        
        //sum of all energy
        MPI_Allreduce(&localEnergy,&globalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        globalEnergy = globalEnergy/2.00;
        if(worldRank==0 && step==1)
            printf("Total energy from algo1 %d = %f\n",worldRank,globalEnergy);
        
        //------------------------------------end of algorithm 1-------------------------------------------------------------------------	
        
        printf("Proc %d: Computational time of algo1: %f secs\n",worldRank,seconds1);
        printf("Proc %d: Computational time of algo2: %f secs\n",worldRank,seconds3);
        MPI_Finalize();
        return 0;
    }
