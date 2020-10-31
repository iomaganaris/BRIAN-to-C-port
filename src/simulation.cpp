/**
 * @file simulation.cpp
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2018
 * @brief Simulation main file.
 *
 * Contains the main function of the program, with all the initializations
 * of the variables needed for a specific simulation.
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "constants.hpp"
#include "synapses.hpp"


//#define NxM
#define MxM
//#define NxMxM	//not working probably
#define allconnected


int main(void){

    std::ofstream spikes, neurons_I, array_A;
    spikes.open("Spikes.txt");
    neurons_I.open("Neurons_I.txt");
    array_A.open("array_A.txt");
    std::ifstream connections;
    connections.open("connections.txt");
	int nruns = 1;

	for(int nrun = 0; nrun < nruns; nrun++){
		double realtime = 0;
		double stime = 0.005; //second
		double stime2 = 50; //second

		double resolution_export = 10 * 1e-3; //every x ms

		int N = 100;
		#ifdef MxM
			N = 0;
		#endif
		int M = 10;

		int N_S;//100;
		int N_Group_S;
		int N_Group_T;

		#ifdef NxM
			N_S = N;//100;
			N_Group_S = 0;
			N_Group_T = M;	//for the simulation we have, normally is N
	    	#endif

		#ifdef MxM
			N_S = 0;//100;
			N_Group_S = M;
			N_Group_T = M;	//for the simulation we have, normally is N
		#endif
	    /*eqs_neuron = """
        dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-vt)/DeltaT)+I-x)/C : volt
        dvt/dt=-(vt-vtrest)/tauvt : volt
        dx/dt=(c*(vm-EL)-x)/tauw : amp #In the standard formulation x is w
        I : amp
    	"""*/

		double input1_pos = 25;
		double input2_pos = 75;
		double rad = 5;
        /* TODO: Handle poisson neurons later
	    //Define input 1
	    double *F_input1;
	    F_input1 = (double*)malloc(sizeof(double)*N_S);
	    //F_input1[input1_pos-rad:input1_pos+rad] = Fon
	    //printf("F_input1\n");
	    for(int i = 0; i < N_S; i++){
	    	F_input1[i] = Foff;			//maybe is not needed
	    	F_input1[i] = exp(-(pow((i+1)-input1_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    	//printf("%lf ",F_input1[i]);
	    }
	    //printf("\n");

	    //Define input 2
	    double *F_input2;
	    F_input2 = (double*)malloc(sizeof(double)*N_S);
	    //F_input2[input2_pos-rad:input2_pos+rad] = Fon
	    for(int i = 0; i < N_S; i++){
	    	F_input2[i] = Foff;			//maybe is not needed
	    	F_input2[i] = exp(-(pow((i+1)-input2_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    }

	    
	    Poisson *input;
	    input = (Poisson*)malloc(sizeof(Poisson)*N_S);
	    for(int i = 0; i<N_S; i++){				// Initialization of Poisson Neurons
	    	input[i].GaussArray = F_input1[i];
	    	input[i].Spike = 0;
	    }
	    */
        std::vector<double> init_vtrest(N_Group_T, vtrest);
        std::vector<double> init_vm(N_Group_T, vtrest+ 0.005); //EL);
        std::vector<double> init_I(N_Group_T, 0);
        std::vector<double> init_x(N_Group_T, 0);
	    AdEx AdEx_neurons(init_vtrest, init_vm, init_I, init_x);

	    std::vector<double> times;
	    std::vector<int> spike_ids;
        int timesteps = (int)(stime/defaultclock_dt);
	    for(auto t = 0; t < timesteps; t++) {
            for (auto i = 0; i < N_Group_S; i++) {
//                if (i % 2 == t % 2) {
//                    times.push_back(t*defaultclock_dt);
//                    spike_ids.push_back(i);
//                }
                times.push_back(t*defaultclock_dt);
                spike_ids.push_back(1);
            }
        }
        Inputs input_neurons(N_Group_S, spike_ids, times);

	    std::vector<double> init_AdEx_FBp(N_Group_T*N_Group_T, 0);
        std::vector<double> init_AdEx_FBn(N_Group_T*N_Group_T, 0);
        std::vector<double> init_AdEx_R(N_Group_T*N_Group_T, 1);
        std::vector<double> init_AdEx_U(N_Group_T*N_Group_T);
        std::vector<double> init_AdEx_A(N_Group_T*N_Group_T);
        int init_const = 0;
        int index = 0;
        for(int i = 0; i < N_Group_T; i++){
            for(int j = 0; j < N_Group_T; j++){
                index = i*N_Group_T + j;
                init_AdEx_U[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin*10;	// takes time
                init_AdEx_A[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin*10;	// takes time
                init_const++;
            }
        }
        Synapses AdEx_AdEx_synapses(AdEx_neurons, AdEx_neurons, init_AdEx_FBp, init_AdEx_FBn, init_AdEx_R, init_AdEx_U, init_AdEx_A);
	    // Initialization of Synapses for external input
        std::vector<double> init_input_FBp(N_Group_S*N_Group_T, 0);
        std::vector<double> init_input_FBn(N_Group_S*N_Group_T, 0);
        std::vector<double> init_input_R(N_Group_S*N_Group_T, 1);
        std::vector<double> init_input_U(N_Group_S*N_Group_T);
        std::vector<double> init_input_A(N_Group_S*N_Group_T);
        for(int i = 0; i < N_Group_S; i++){
            for(int j = 0; j < N_Group_T; j++){
                index = i*N_Group_T + j;
                init_input_U[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
                init_input_A[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
                init_const++;
            }
        }
        Synapses Input_AdEx_synapses(input_neurons, AdEx_neurons, init_input_FBp, init_input_FBn, init_input_R, init_input_U, init_input_A);
		std::cout << "timesteps=" << timesteps << std::endl;
		for(int t = 0; t < timesteps; t++){			//add monitors for the variables we care about
			input_neurons.generate_spikes(t*defaultclock_dt);
			AdEx_neurons.solve_neurons();
            //input_neurons.print_spikes();
            AdEx_neurons.print_spikes();
            AdEx_neurons.print_neurons();
			//Input_AdEx_synapses.UpdateSynapses_pre(t*defaultclock_dt);
			AdEx_AdEx_synapses.UpdateSynapses_pre(t*defaultclock_dt);
            //AdEx_AdEx_synapses.print_synapses();
            //Input_AdEx_synapses.UpdateSynapses_post(t*defaultclock_dt);
            AdEx_AdEx_synapses.UpdateSynapses_post(t*defaultclock_dt);
		}
	}
    spikes.close();
    neurons_I.close();
    array_A.close();
    connections.close();
	return 0;
}
