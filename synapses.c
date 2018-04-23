/**
 * @file synapses.c
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2014
 * @brief File containing all the functions of Synapses.
 *
 * It contains the function for updating the state of a synapse due to
 * a pre- and post-synaptic event and also a function for printing all the 
 * elements of the 2D matrix of synapses.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "synapses.h"

extern double defaultclock_dt;

//double taum_S = 10 * 1e-3; //ms
extern double Ee;// = 0 * 1e-3; //mV
extern double tuae;// = 2 * 1e-3; //ms
extern double Fon;// = 50; //Hz
extern double Foff;// = 3; //Hz

extern double s;// = 100*1e-10; 
extern double Amax;// = 2.0;
extern double Amin;// = 0;
extern double Ainit;// = 0.1;
extern double Umax;// = 1.0;
extern double Umin;// = 0;
extern double Uinit;// = 0.1;

extern double dFBn;// = 0;
extern double dFBp;// = 0;
extern double dFFp;// = 0;

//#Short-term plasticity params
extern double tau_u;// = 50 * 1e-3;	//ms
extern double tau_r;// = 200 * 1e-3;	//ms

//#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
//double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
extern double AFBn;// = 0.1771;
extern double tau_FBn;// = 0.0327 * 1e3 * 1e-3;	//ms
extern double AFBp;// = 0.1548;
extern double tau_FBp;// = 0.2302 * 1e3 * 1e-3;	//ms
extern double AFFp;// = 0.0618;
extern double tau_FFp;// = 0.0666 * 1e3 * 1e-3;	//ms
//#etaU = 0.35
extern double etaU;// = 0.15;
extern double etaA;// = 0.15;
//#etaA = 0.35    

void UpdateSynapses_pre(Synapse** Synapses, Neuron* neurons, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t){
	//an baloume ekswteriki eisodo thelei i<N+1 logika
	for (int i = 0; i < N_S + N_Group_S; i++){
		if (SpikeArray[i+N_Group_T-N_Group_S] > 0){
			printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
			for (int j = 0; j < N_Group_T; j++){
				//auto borei na ginei pio apodotika alla diskoleuei, gi arxi trwme xwro me ta struct synapses
				if (Synapses[i][j].conn){
					printf("i = %d, j = %d\n",i,j);
					Synapses[i][j].FFp = Synapses[i][j].FFp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FFp);
					Synapses[i][j].FBn = Synapses[i][j].FBn * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBn);
					Synapses[i][j].u = Synapses[i][j].U + (-Synapses[i][j].U + Synapses[i][j].u) * exp(-(-Synapses[i][j].lastupdate + t)/tau_u);
					Synapses[i][j].FBp = Synapses[i][j].FBp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBp);
					Synapses[i][j].R = (Synapses[i][j].R - 1) * exp(-(-Synapses[i][j].lastupdate + t)/tau_r) + 1;
					Synapses[i][j].target_I = s * Synapses[i][j].A * Synapses[i][j].R * Synapses[i][j].u;
					Synapses[i][j].U = Synapses[i][j].U + etaU * (-AFBn * Synapses[i][j].FBn * Synapses[i][j].FBp + AFBp * Synapses[i][j].FBp * Synapses[i][j].FFp);
					if (Synapses[i][j].U < Umin) Synapses[i][j].U = Umin;
					else if (Synapses[i][j].U > Umax) Synapses[i][j].U = Umax;
					Synapses[i][j].w = Synapses[i][j].U * Synapses[i][j].A;
					Synapses[i][j].FFp += 1;
					Synapses[i][j].R -= Synapses[i][j].R * Synapses[i][j].u;
					Synapses[i][j].u += Synapses[i][j].U * (1 - Synapses[i][j].u);
					Synapses[i][j].lastupdate = t;
				}
			}
		}
	}
	for (int j = 0; j < N_Group_T; j++){
		for (int i = 0; i < N_S + N_Group_S; i++){
			if (Synapses[i][j].conn && SpikeArray[i+N_Group_T-N_Group_S]){
				//printf("i: %d, j: %d\n",i,j);
				neurons[j].I = Synapses[i][j].target_I;
				printf("j: %d, neurons[j]= %lf",j,neurons[j].I);
				//break;
			}
		}
	}
}

void UpdateSynapses_post(Synapse** Synapses, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t){
	
	for (int i = 0; i < N_Group_T; i++){
	    if (SpikeArray[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_S + N_Group_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
	            if (Synapses[j][i].conn){
	            	printf("i = %d, j = %d\n",i,j);
	            	fflush(stdout);
	                Synapses[j][i].FFp = Synapses[j][i].FFp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FFp);
	                Synapses[j][i].FBn = Synapses[j][i].FBn * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBn);
	                Synapses[j][i].u = Synapses[j][i].U + (-Synapses[j][i].U + Synapses[j][i].u) * exp(-(-Synapses[j][i].lastupdate + t)/tau_u);
	                Synapses[j][i].FBp = Synapses[j][i].FBp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBp);
	                Synapses[j][i].R = (Synapses[j][i].R - 1) * exp(-(-Synapses[j][i].lastupdate + t)/tau_r) + 1;
	                /*Synapses[j][i].A = Synapses[j][i].A + etaA * (AFFp * Synapses[j][i].FFp * Synapses[j][i].FBn);
	                double mean = 0;
	                int num = 0;
	                for (int k=0; k<N_S; k++){
	                    for (int l=0; l<N_T; l++){
	                        //auto edw to if boroume na to apofigoume an arxikopoioume tis metablites tou struct sto 0
	                        //if (Synapses[k][l].conn && (SpikeArray[k] || SpikeArray[l])){
	                    	if (Synapses[k][l].conn){
	                            mean += AFFp * Synapses[k][l].FFp * Synapses[k][l].FBn;
	                            num++;
	                        }
	                    }
	                }
	                mean = (double)mean / num;
	                printf("mean = %lf",mean);
	                Synapses[j][i].A = Synapses[j][i].A - etaA * 0.5 * mean; //amfibola swsto, sigoura mi apodotiko
	                if (Synapses[j][i].A < Amin) Synapses[j][i].A = Amin;
	                else if (Synapses[j][i].A > Amax) Synapses[j][i].A = Amax;
	                Synapses[j][i].w = Synapses[j][i].U * Synapses[j][i].A;
	                Synapses[j][i].FBp += 1;
	                Synapses[j][i].FBn += 1;*/

	            }
	        }
	    }
	}
	for (int i = 0; i < N_Group_T; i++){
	    if (SpikeArray[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_S + N_Group_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
	            if (Synapses[j][i].conn){
	            	Synapses[j][i].A = Synapses[j][i].A + etaA * (AFFp * Synapses[j][i].FFp * Synapses[j][i].FBn);
	                long double mean = 0;
	                int num = 0;
	                for (int k=0; k<N_S+N_Group_S; k++){
	                    for (int l=0; l<N_Group_T; l++){
	                        //auto edw to if boroume na to apofigoume an arxikopoioume tis metablites tou struct sto 0
	                        //if (Synapses[k][l].conn && (SpikeArray[k] || SpikeArray[l])){
	                    	if (Synapses[k][l].conn && SpikeArray[l]){
	                            mean += AFFp * Synapses[k][l].FFp * Synapses[k][l].FBn;
	                            num++;
	                        }
	                    }
	                }
	                mean = (double)mean / (double)num;
	                printf("mean = %Le",mean);
	                Synapses[j][i].A = Synapses[j][i].A - etaA * 0.5 * mean; //amfibola swsto, sigoura mi apodotiko
	                if (Synapses[j][i].A < Amin) Synapses[j][i].A = Amin;
	                else if (Synapses[j][i].A > Amax) Synapses[j][i].A = Amax;
	            }
	        }
	    }
	}
	for (int i = 0; i < N_Group_T; i++){
	    if (SpikeArray[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_S + N_Group_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
	            if (Synapses[j][i].conn){
	            	Synapses[j][i].w = Synapses[j][i].U * Synapses[j][i].A;
	                Synapses[j][i].FBp += 1;
	                Synapses[j][i].FBn += 1;
	                Synapses[j][i].lastupdate = t;
	            }
	        }
	    }
	}


}

void print_synapses(Synapse** syn, int N_S, int N_T){
	printf("conn\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			//printf("conn= %d, w= %lf, FFp= %lf, FBp= %lf, FBn= %lf, R= %lf, u= %lf, U= %lf, A= %lf, lastup= %lf, target_I= %lf\n",syn[i][j].conn,syn[i][j].w,syn[i][j].FFp,syn[i][j].FBp,syn[i][j].FBn,syn[i][j].R,syn[i][j].u,syn[i][j].U,syn[i][j].A,syn[i][j].lastupdate,syn[i][j].target_I);
			printf("%d ", syn[i][j].conn);
			//if((i*N_T+j)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nw\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].w);
			if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		//printf("\n");		
	}
	printf("\nFFp\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FFp);
		printf("\n");	
	}
	printf("\nFBp\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FBp);	
		printf("\n");
	}
	printf("\nFBn\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FBn);
		printf("\n");		
	}
	printf("\nR\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].R);
		printf("\n");
	}
	printf("\nu\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].u);
			if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		//printf("\n");		
	}
	printf("\nU\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].U);
			if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		//printf("\n");
	}
	printf("\nA\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].A);	
			if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		//printf("\n");
	}
	printf("\nlastupdate\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%lf, ", syn[i][j].lastupdate);		
		//printf("\n");
	}
	printf("\ntarget_I\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].target_I);
		}
		printf("\n");	
	}
	printf("\n");
}