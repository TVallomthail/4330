/******************************************************************************
* Final project program
* Group Members: Treasa Vallomthail, Andrew Becker, Michael Toelle
* Last revised: 12/1/2022
* DESCRIPTION: Serial Simulation of the SEIR model using 2D array processing in
*                conjunction with MPI.
******************************************************************************/
#include <iostream>
#include <math.h>
#include <string>
#include <time.h>

using namespace std;

int x, y;
double Scount = 480.0;
double Icount = 326.0;
double Rcount = 194.0;
int Population[100][100];

void initPopGrid(int pop) {
    double popRoot = sqrt(pop);
    x = ceil(popRoot);
    y = ceil(popRoot);

    //initialize people
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            int randNum = (rand() % (3 - 1 + 1)) + 1;
            Population[i][j] = randNum;
        }
    }
}

void doWork(int x, int y) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            if (Population[i][j] == 1) { //1 stands for Susceptible
                Scount++;
            }
            else if (Population[i][j] == 2) { //2 stands for Infected
                Icount++;
            }
            else {
                Rcount++; //else Recovered
            }
        }
    }
}

/*Calculations adapted from Brian Sullivan's implementation for
SIR Model for Disease Spread */
void calculate(int days, int population) {

    double dt = 1.0; //time step in days
    double beta = 0.2; //infection rate
    double gamma = 0.07; //recovery rate

    double S[480], I[326], R[194];

    I[0] = Icount / population; //initial infective population
    S[0] = Scount / population - I[0]; //initial susceptible population
    R[0] = Rcount / population; //initial recovered population

    //print initials
    cout << "Fraction of population susceptible to infection at beginning of observation: " << S[0] << endl;
    cout << "Fraction of population infected at beginning of observation: " << I[0] << endl;
    cout << "Fraction of population recovered at beginning of observation: " << R[0] << endl;
    

    for (int i = 0; i < days; i++) {
        S[i + 1] = S[i] - beta * (S[i] * I[i]) * dt;
        I[i + 1] = I[i] + (beta * S[i] * I[i] - gamma * I[i]) * dt;
        R[i + 1] = R[i] + gamma * I[i] * dt;

        //print values
        cout << "Fraction of population susceptible to infection at day " << i << ": " << S[i] << endl;
        cout << "Fraction of population infected at day " << i << ": " << I[i] << endl;
        cout << "Fraction of population recovered at day " << i << ": " << R[i] << endl;
    }
}


//create master thread and worker threads to find population percentages
//do a loop of calculations afterward to simulate a week (with a simple loop and the equations)

int main(int argc, char* argv[])
{
    int population = 1000;
    int days = 30;
    cout << "Using population size of 1000" << endl;
    cout << "Over 30 days" << endl;

    initPopGrid(population);

    srand(time(0));
    doWork(x, y);
    calculate(days, population);

    return 0;
}