#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "neuron.h"

extern double defaultclock_dt;// = 1*1e-3;	//ms

//# Adex Parameters
extern double C;// = 281*1e-12;	//pF
extern double gL;// = 30*1e-9; //nS
extern double taum;// = 281*1e-12 / 30*1e-9;	// C/gL	// double initilization of taum(?)
extern double EL;// = -70.6*1e-3;	//mV
extern double DeltaT;// = 2*1e-3;	//mV
extern double vti;// = -50.4*1e-3;	//mV
//#vtrest = vti + 5 * DeltaT
extern double vtrest;// = -45*1e-3;	//mV
extern double VTmax;// = 18*1e-3;	//mV
extern double tauvt;// = 50*1e-3;	//ms

extern double tauw;// = 144*1e-3;	//ms
extern double c;// = 4*1e-9;	//ns
extern double b;// = 0.0805*1e-9;	//nA
extern double Vr;// = -70.6*1e-3;	//mV

void print_neurons(Neuron* neurons, int N){
	printf("vt\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].vt);
	}
	printf("\nvm\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].vm);
	}
	printf("\nI\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].I);
	}
	printf("\nx\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].x);
	}
	printf("\nSpike\n");
	for(int i = 0; i < N; i++){
		printf("%d, ",neurons[i].Spike);
	}
	printf("\n");
}

void resetNeuron(Neuron* neuron) {
    neuron->vm = Vr;
    neuron->x += b;
    neuron->vt = VTmax;
}

double random_0_1(){
	return (double)rand()/(double)((unsigned)RAND_MAX+1);
}

void PoissonThreshold(Poisson* input, int N_S, int N_T, int* SpikeArray){
	for(int i = 0; i<N_S; i++){
		if(random_0_1()<input[i].GaussArray*defaultclock_dt){
			input[i].Spike = 1;
			SpikeArray[N_T+i] = 1;
			printf("PoissonThreshold i= %d\n",i);
		} 
		else{
			input[i].Spike = 0;
			SpikeArray[N_T+i] = 0;
		}
	}
	return;
}

void SolveNeurons(Neuron* neurons, int N, int *SpikeArray){
    for(int i = 0; i < N; i++){
        /*if (rand() % 10 < 5) {
            SpikeArray[i] = 0;
        }
        else {
            SpikeArray[i] = 1;
        }*/
        double _vm, _vt, _x;
        _vm = (gL*(EL-neurons[i].vm)+gL*DeltaT*exp((neurons[i].vm-neurons[i].vt)/DeltaT)+neurons[i].I-neurons[i].x)/C;
        //printf("_vm = %.20f\n", _vm);
        _vt = -(neurons[i].vt-vtrest)/tauvt;
        //printf("_vt = %.20f\n", _vm);
        _x = (c*(neurons[i].vm-EL)-neurons[i].x)/tauw;
        //printf("_x = %.20f\n", _x);
        neurons[i].vm += _vm * defaultclock_dt;
        neurons[i].vt += _vt * defaultclock_dt;
        neurons[i].x += _x * defaultclock_dt;
        if(neurons[i].vm > neurons[i].vt){
            //printf("Reset\n");
            resetNeuron(&neurons[i]);
            SpikeArray[i] = 1;
        }
        else SpikeArray[i] = 0;
    }
}