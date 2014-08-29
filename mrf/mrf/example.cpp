#include "mrf.h"
#include "ICM.h"
#include "GCoptimization.h"
#include "MaxProdBP.h"


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>


const int sizeX = 200;
const int sizeY = 200;
const int K = 16;

MRF::CostVal D[sizeX*sizeY*K];
MRF::CostVal V[K*K];
MRF::CostVal hCue[sizeX*sizeY];
MRF::CostVal vCue[sizeX*sizeY];

EnergyFunction* generate_DataARRAY_SmoothFIXED_FUNCTION()
{
    int i, j;


    // generate function
    for (i=0; i<K; i++)
    for (j=i; j<K; j++)
    {
        V[i*K+j] = V[j*K+i] = (i == j) ? 0 : 2;
    }
    MRF::CostVal* ptr;
    for (ptr=&D[0]; ptr<&D[sizeX*sizeY*K]; ptr++) *ptr = rand() % 10;
    for (ptr=&hCue[0]; ptr<&hCue[sizeX*sizeY]; ptr++) *ptr = rand() % 3 - 1; // negative multiplier possible
    for (ptr=&vCue[0]; ptr<&vCue[sizeX*sizeY]; ptr++) *ptr = rand() % 3;

    // allocate eng
    DataCost *data         = new DataCost(D);
    SmoothnessCost *smooth = new SmoothnessCost(V,hCue,vCue);
    EnergyFunction *eng    = new EnergyFunction(data,smooth);

    return eng;
}

EnergyFunction* generate_DataARRAY_SmoothTRUNCATED_LINEAR()
{
    // generate function
    MRF::CostVal* ptr;
    for (ptr=&D[0]; ptr<&D[sizeX*sizeY*K]; ptr++) *ptr = rand() % 10;
    for (ptr=&hCue[0]; ptr<&hCue[sizeX*sizeY]; ptr++) *ptr = rand() % 3;
    for (ptr=&vCue[0]; ptr<&vCue[sizeX*sizeY]; ptr++) *ptr = rand() % 3;
    MRF::CostVal smoothMax = 5, lambda = 2;

    // allocate eng

    DataCost *data         = new DataCost(D);
    SmoothnessCost *smooth = new SmoothnessCost(1,smoothMax,lambda,hCue,vCue);
    EnergyFunction *eng    = new EnergyFunction(data,smooth);

    return eng;
}


EnergyFunction* generate_DataARRAY_SmoothTRUNCATED_QUADRATIC()
{
    
    // generate function
    MRF::CostVal* ptr;
    for (ptr=&D[0]; ptr<&D[sizeX*sizeY*K]; ptr++) *ptr = rand() % 10;
    for (ptr=&hCue[0]; ptr<&hCue[sizeX*sizeY]; ptr++) *ptr = rand() % 3;
    for (ptr=&vCue[0]; ptr<&vCue[sizeX*sizeY]; ptr++) *ptr = rand() % 3;
    MRF::CostVal smoothMax = 5,lambda=2;

    // allocate eng
    DataCost *data         = new DataCost(D);
    SmoothnessCost *smooth = new SmoothnessCost(2,smoothMax,lambda,hCue,vCue);
    EnergyFunction *eng    = new EnergyFunction(data,smooth);

    return eng;
}


MRF::CostVal dCost(int pix, int i)
{
    return (pix*i + i + pix) % 10;
}

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{
    return (pix1*(i+1)*(j+2) + pix2*i*j*pix1 - 2*i*j*pix1) % 10;
}


EnergyFunction* generate_DataFUNCTION_SmoothGENERAL_FUNCTION()
{
    DataCost *data         = new DataCost(dCost);
    SmoothnessCost *smooth = new SmoothnessCost(fnCost);
    EnergyFunction *eng    = new EnergyFunction(data,smooth);

    return eng;
}

int main()
{
    MRF* mrf;
    EnergyFunction *eng;
    MRF::EnergyVal E;
    float t,tot_t;
    int iter;

    int seed = 1124285485;
    srand(seed);

    // There are 4 sample energies below to play with. Uncomment 1 at a time 

    //eng = generate_DataARRAY_SmoothFIXED_FUNCTION();
    //eng = generate_DataARRAY_SmoothTRUNCATED_LINEAR();
    eng = generate_DataARRAY_SmoothTRUNCATED_QUADRATIC();
    //eng = generate_DataFUNCTION_SmoothGENERAL_FUNCTION();



    ////////////////////////////////////////////////
    //                     ICM                    //
    ////////////////////////////////////////////////
    printf("\n*******Started ICM *****\n");

    mrf = new ICM(sizeX,sizeY,K,eng);
    mrf->initialize();
    mrf->clearAnswer();
    
    E = mrf->totalEnergy();
    printf("Energy at the Start= %d (%d,%d)\n", E,mrf->smoothnessEnergy(),mrf->dataEnergy());

    tot_t = 0;
    for (iter=0; iter<6; iter++)
    {
        mrf->optimize(10, t);

        E = mrf->totalEnergy();
        tot_t = tot_t + t ;
        printf("energy = %d (%f secs)\n", E, tot_t);
    }

    delete mrf;

    ////////////////////////////////////////////////
    //          Graph-cuts expansion              //
    ////////////////////////////////////////////////
    printf("\n*******Started the graph-cuts expansion *****\n");
    mrf = new Expansion(sizeX,sizeY,K,eng);
    mrf->initialize();
    mrf->clearAnswer();
    
    E = mrf->totalEnergy();
    printf("Energy at the Start= %d (%d,%d)\n", E,mrf->smoothnessEnergy(),mrf->dataEnergy());

    tot_t = 0;
    for (iter=0; iter<6; iter++)
    {
        mrf->optimize(1, t);

        E = mrf->totalEnergy();
        tot_t = tot_t + t ;
        printf("energy = %d (%f secs)\n", E, tot_t);
    }

    delete mrf;

    ////////////////////////////////////////////////
    //          Graph-cuts swap                   //
    ////////////////////////////////////////////////

    printf("\n*******Started the graph-cuts swap *****\n");
    mrf = new Swap(sizeX,sizeY,K,eng);
    mrf->initialize();
    mrf->clearAnswer();
    
    E = mrf->totalEnergy();
    printf("Energy at the Start= %d (%d,%d)\n", E,mrf->smoothnessEnergy(),mrf->dataEnergy());

    tot_t = 0;
    for (iter=0; iter<6; iter++)
    {
        mrf->optimize(1, t);

        E = mrf->totalEnergy();
        tot_t = tot_t + t ;
        printf("energy = %d (%f secs)\n", E, tot_t);
    }

    
    delete mrf;


    ////////////////////////////////////////////////
    //          Belief Propagation                //
    ////////////////////////////////////////////////

    printf("\n*******  Started MaxProd Belief Propagation *****\n");
    mrf = new MaxProdBP(sizeX,sizeY,K,eng);
    mrf->initialize();
    mrf->clearAnswer();
    
    E = mrf->totalEnergy();
    printf("Energy at the Start= %d (%d,%d)\n", E,mrf->smoothnessEnergy(),mrf->dataEnergy());

    tot_t = 0;
    for (iter=0; iter < 10; iter++)
    {
        mrf->optimize(1, t);

        E = mrf->totalEnergy();
        tot_t = tot_t + t ;
        printf("energy = %d (%f secs)\n", E, tot_t);
    }

    
    delete mrf;

    return 0;
}

