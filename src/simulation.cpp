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

#include "constants.hpp"
#include "synapses.hpp"

int main(int argc, char **argv){

    if (argc < 4)
    {
        std::cout << "Simulation needs 3 arguments. Simulation time, number of Input neurons and number of AdEx neurons.\n";
        return 1;
    }

    double stime = std::stof(argv[1]); // seconds

    int N_Group_S = std::stoi(argv[2]);
    int N_Group_T = std::stoi(argv[3]);

    std::cout << "Simulating " << N_Group_S << " Input and " << N_Group_T << " AdEx Neurons\n";

    double input1_pos = 25;
    double rad = 5;

    /// Initializing AdEx Neurons
    std::vector<double> init_vtrest(N_Group_T, vtrest);
    std::vector<double> init_vm(N_Group_T, vtrest+ 0.005);
    std::vector<double> init_I(N_Group_T, 0);
    std::vector<double> init_x(N_Group_T, 0);
    AdEx AdEx_neurons(init_vtrest, init_vm, init_I, init_x);

    /// Initializing Input Neurons firing every 2 time steps
    std::vector<double> times;
    std::vector<int> spike_ids;
    int timesteps = (int)(stime/defaultclock_dt);
    for(auto t = 0; t < timesteps; t++) {
        for (auto i = 0; i < N_Group_S; i++) {
            if (t % 2 == 0) {
                times.push_back(t*defaultclock_dt);
                spike_ids.push_back(i);
            }
        }
    }
    Inputs input_neurons(N_Group_S, spike_ids, times);

    /// Initializing synapses between AdEx Neurons
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

    /// Initializing synapses between Input and AdEx Neurons
    std::vector<double> init_input_FBp(N_Group_S*N_Group_T, 0);
    std::vector<double> init_input_FBn(N_Group_S*N_Group_T, 0);
    std::vector<double> init_input_R(N_Group_S*N_Group_T, 1);
    std::vector<double> init_input_U(N_Group_S*N_Group_T);
    std::vector<double> init_input_A(N_Group_S*N_Group_T);
    init_const = 0;
    for(int i = 0; i < N_Group_S; i++){
        for(int j = 0; j < N_Group_T; j++){
            index = i*N_Group_T + j;
            init_input_U[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
            init_input_A[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
            init_const++;
        }
    }
    Synapses Input_AdEx_synapses(input_neurons, AdEx_neurons, init_input_FBp, init_input_FBn, init_input_R, init_input_U, init_input_A);

    /// Main simulation loop
    std::cout << "timesteps=" << timesteps << "\n";
    for(int t = 0; t < timesteps; t++){
        input_neurons.generate_spikes(t*defaultclock_dt);
        AdEx_neurons.solve_neurons();
        //input_neurons.print_spikes();
        AdEx_neurons.print_spikes();
        //AdEx_neurons.print_neurons();
        Input_AdEx_synapses.UpdateSynapses_pre(t*defaultclock_dt);
        AdEx_AdEx_synapses.UpdateSynapses_pre(t*defaultclock_dt);
        //AdEx_AdEx_synapses.print_synapses();
        Input_AdEx_synapses.UpdateSynapses_post(t*defaultclock_dt);
        AdEx_AdEx_synapses.UpdateSynapses_post(t*defaultclock_dt);
        //Input_AdEx_synapses.print_synapses();
    }

	return 0;
}
