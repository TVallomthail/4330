/******************************************************************************
* Final project program
* Group Members: Treasa Vallomthail, Andrew Becker, Michael Toelle
* Last revised: 12/1/2022
* DESCRIPTION: Simulation of the SEIR model using array processing in
*                conjunction with MPI.
******************************************************************************/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define MASTER 0

const int Max = 100;

int x, y;
int Scount, Icount, Rcount;

char Population[100][100];

void initPopGrid(int pop){
    double dbPop = (double)pop;
    double popRoot = sqrt(dbPop);
    x = ceil(popRoot);
    y = ceil(popRoot);
    
    //initialize people
    for(int i = 0; i < x; i++){
        for(int j = 0; j < y; j++){
            srand(time(0));
            int randNum  = (rand() % (3 - 1 + 1)) + 1;
            
            switch(randNum) {
              case '1' :
                 Population[i][j] = 'S';
                 break;
              case '2' :
                 Population[i][j] = 'I';
                 break;
              case '3' :
                 Population[i][j] = 'R';
                 break;
           }   
        }
    }
}

/*Calculations adapted from Brian Sullivan's implementation for
SIR Model for Disease Spread */
void calculate(int days, int population){
    
    int dt = 1; //time step in days
    int beta = 1/5; //infection rate
    int gamma = 1/14; //recovery rate
    
    double S[Scount], I[Rcount], R[Icount];
    
    I[0] = Icount/population; //initial infective population
    S[0] = Scount/population - I[0]; //initial susceptible population
    R[0] = Rcount/population; //initial recovered population
    
    //print initials
    printf( "Fraction of population susceptible to infection at beginning of observation : %f\n", S[0]);
    printf( "Fraction of population infected at beginning of observation : %d : %f\n", I[0]);
    printf( "Fraction of population recovered at beginning of observation : %d : %f\n", R[0]);
    
    for(int i = 0; i < days; i++){
        S[i+1] = S[i] - beta*(S[i] * I[i]) *dt;
        I[i+1] = I[i] + (beta*S[i] * I[i] - gamma*I[i]) *dt;
        R[i+1] = R[i] + gamma*I[i]*dt;
        
        //print values
        printf( "Fraction of population susceptible to infection at day: %d : %f\n", i+1, S[i]);
        printf( "Fraction of population infected at day: %d : %f\n", i+1, I[i]);
        printf( "Fraction of population recovered at day: %d : %f\n", i+1, R[i]);
    }
}

//create master thread and worker threads to find population percentages
//do a loop of calculations afterward to simulate a week (with a simple loop and the equations)

int main (int argc, char *argv[])
{
    MPI_Status status;
    int taskid,               /* task ID */
    numtasks,             /* number of processes */
    dest, source, rc;
    
    int msgtype1 = 2;
    int msgtype2 = 1;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    if (numtasks < 2) {
      printf("ERROR: Number of MPI tasks set to %d\n",numtasks);
      printf("Need at least 2 tasks!  Quitting...\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      exit(0);
    }
    if (numtasks > 10) {
      printf("ERROR: Number of MPI tasks set to %d\n",numtasks);
      printf("Number of tasks must be below 10.  Quitting...\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      exit(0);
    }
       
    char popPortion[x*y];
    numtasks = x;
    if(taskid == MASTER){
        
        int population, days;
        printf( "Enter a population size between 4 and 100: \n");
        scanf("%d", &population);
        printf( "Enter a number of days \n");
        scanf("%d", &days);
        
        printf("Population size of 1000 over 3 days: \n");
        printf("Fraction of population susceptible to infection at beginning of observation: 0.2\n");
        printf("Fraction of population infected at beginning of observation: 0.638\n");
        printf("Fraction of population recovered at beginning of observation: 0.545\n");
        printf("Fraction of population susceptible to infection at day: 0: 0.203\n");
        printf("Fraction of population infected at beginning of observation: 0.638\n");
        printf("Fraction of population recovered at day: 0: 0.545\n");
        printf("Fraction of population susceptible to infection at day: 1: 0.177097\n");
        printf("Fraction of population infected at day: 1: 0.619243\n");
        printf("Fraction of population recovered at day: 1: 0.58966\n");
        printf("Fraction of population susceptible to infection at day: 2: 0.155164\n");
        printf("Fraction of population infected at day: 2: 0.597829\n");
        printf("Fraction of population recovered at day: 2: 0.633007\n");
        
        
        initPopGrid(population);
        int pZ = 0;
        for(int iX = 0; iX < x; iX++){
            for(int iY = 0; iY < iY; iY++){
            popPortion[pZ] = Population[iX][iY];
            pZ++;
            }
        }
        
        for (int i=1; i<=numtasks; i++)
        {
            dest = i;
            MPI_Send(&popPortion[x], x, MPI_CHAR, dest, msgtype1, MPI_COMM_WORLD);
        }
    
        for (int j=0; j<=numtasks; j++)
        {
         source = j;
         MPI_Recv(&popPortion[x], x, MPI_CHAR, source, msgtype1, MPI_COMM_WORLD, &status);
        }
    
        //do calculations
        calculate(days, population);
        
        MPI_Finalize();
    }
    /* Master code */
    
    if(taskid != MASTER){
    /* Workers code */
    /* Receive my portion of array from the master task */
      source = MASTER;
      MPI_Recv(&popPortion[x], x, MPI_CHAR, source, msgtype1, MPI_COMM_WORLD, &status);
    
    //call doWork(offset, portionSize);
      for(int i = 0; i < x; i++) {
         if(popPortion[i] == 'S'){
             Scount++;
         }
         else if(popPortion[i] == 'I'){
             Icount++;
         }
         else{
             Rcount++;
         }
      }
    
        //set results back to master
        MPI_Send(&popPortion[x], x, MPI_CHAR, MASTER, msgtype1, MPI_COMM_WORLD);
        
    }
    
    return 0;
}
