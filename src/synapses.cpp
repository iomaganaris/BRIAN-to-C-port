/**
 * @file synapses.cpp
 * @author Ioannis Magkanaris
 * @date 4 November 2020
 * @brief File containing all the functions of Synapses
 *
 * It contains the function for updating the state of a synapse due to
 * a pre- and post-synaptic event and also a function for printing all the 
 * elements of the 2D matrix of synapses
 */

#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "synapses.hpp"

void Synapses::init() {
    w.resize(N_Group_S*N_Group_T, 0);
    FFp.resize(N_Group_S*N_Group_T, 0);
    u.resize(N_Group_S*N_Group_T, 1);
    lastupdate.resize(N_Group_S*N_Group_T, 0);
    target_I.resize(N_Group_S*N_Group_T, 0);
}

void Synapses::UpdateSynapses_pre(double t){
    auto pre_spikes = pre_neurons.get_spikes();
    auto index = 0;
	for (int i = 0; i < N_Group_S; i++){
		if (pre_spikes[i] > 0){
			//printf("i= %d SpikeArray[i]= %d\n",i,pre_spikes[i]);
			for (int j = 0; j < N_Group_T; j++){
                //printf("i = %d, j = %d\n",i,j);
                index = i*N_Group_T+j;
                FFp[index] = FFp[index] * exp(-(-lastupdate[index] + t)/tau_FFp);
                FBn[index] = FBn[index] * exp(-(-lastupdate[index] + t)/tau_FBn);
                u[index] = U[index] + (-U[index] + u[index]) * exp(-(-lastupdate[index] + t)/tau_u);
                FBp[index] = FBp[index] * exp(-(-lastupdate[index] + t)/tau_FBp);
                R[index] = (R[index] - 1) * exp(-(-lastupdate[index] + t)/tau_r) + 1;
                target_I[index] = s * A[index] * R[index] * u[index];
                U[index] = U[index] + etaU * (-AFBn * FBn[index] * FBp[index] + AFBp * FBp[index] * FFp[index]);
                if (U[index] < Umin) U[index] = Umin;
                else if (U[index] > Umax) U[index] = Umax;
                w[index] = U[index] * A[index];
                FFp[index] += 1;
                R[index] -= R[index] * u[index];
                u[index] += U[index] * (1 - u[index]);
                lastupdate[index] = t;
			}
		}
	}

	post_neurons.update_I(pre_spikes, target_I);
	/* This should be done in update_I
    for (int j = 0; j < N_Group_T; j++){
		for (int i = 0; i < N_Group_S; i++){
			if (pre_spikes[i]){
				//printf("i: %d, j: %d\n",i,j);
                post_neurons->I[j] = target_I[i*N_Group_T+j];
				printf("j: %d, neurons[j]= %lf",j,post_neurons->I[j]);
				//break;
			}
		}
	}*/
}

void Synapses::UpdateSynapses_post(double t){
	auto post_spikes = post_neurons.get_spikes();
    auto index = 0;
	for (int i = 0; i < N_Group_T; i++){
	    if (post_spikes[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_Group_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
                //printf("i = %d, j = %d\n",i,j);
                //fflush(stdout);
                index = i+j*N_Group_T;
                FFp[index] = FFp[index] * exp(-(-lastupdate[index] + t)/tau_FFp);
                FBn[index] = FBn[index] * exp(-(-lastupdate[index] + t)/tau_FBn);
                u[index] = U[index] + (-U[index] + u[index]) * exp(-(-lastupdate[index] + t)/tau_u);
                FBp[index] = FBp[index] * exp(-(-lastupdate[index] + t)/tau_FBp);
                R[index] = (R[index] - 1) * exp(-(-lastupdate[index] + t)/tau_r) + 1;
                A[index] = A[index] + etaA * (AFFp * FFp[index] * FBn[index]);
                double mean = 0;
                for (int k=0; k<N_Group_S; k++){
                    for (int l=0; l<N_Group_T; l++){
                        //auto edw to if boroume na to apofigoume an arxikopoioume tis metablites tou struct sto 0
                        //if (Synapses[k][l].conn && (SpikeArray[k] || SpikeArray[l])){
                        mean += AFFp * FFp[index] * FBn[index];
                    }
                }
                mean = (double)mean / (N_Group_S*N_Group_T);
                //printf("mean = %lf",mean);
                A[index] = A[index] - etaA * 0.5 * mean; //amfibola swsto, sigoura mi apodotiko
                if (A[index] < Amin) A[index] = Amin;
                else if (A[index] > Amax) A[index] = Amax;
                w[index] = U[index] * A[index];
                FBp[index] += 1;
                FBn[index] += 1;
	        }
	    }
	}
	/*
	for (int i = 0; i < N_Group_T; i++){
	    if (post_neurons.spikes[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
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
	}*/
}

void Synapses::print_synapses(){
    std::cout << "\nw\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
			std::cout << w[i*N_Group_T+j] << ", ";
		std::cout << "\n";
	}
    std::cout << "\nFFp\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << FFp[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nFBp\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << FBp[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nFBn\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << FBn[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nR\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << R[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nu\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << u[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nU\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << U[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nA\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << A[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\nlastupdate\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << lastupdate[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\ntarget_I\n";
	for(int i =0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++)
            std::cout << target_I[i*N_Group_T+j] << ", ";
        std::cout << "\n";
	}
    std::cout << "\n";
}