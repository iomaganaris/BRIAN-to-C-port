/**
 * @file simulation.cpp
 * @author Ioannis Magkanaris
 * @date 4 November 2020
 * @brief Simulation main file.
 *
 * Contains the initialization of the circuit and the main simulation loop.
 */
#include <iostream>
#include <random>

#include "constants.hpp"
#include "synapses.hpp"

int main(int argc, char **argv){

    if (argc < 4)
    {
        std::cout << "Simulation needs 3 arguments. Simulation time (s), number of Input neurons and number of AdEx neurons.\n";
        return 1;
    }

    /// Simulation time in seconds
    const double stime = std::stof(argv[1]);
    /// Number of input neurons
    const int N_Group_S = std::stoi(argv[2]);
    /// Number of AdEx neurons
    const int N_Group_T = std::stoi(argv[3]);

    std::cout << "Simulating " << N_Group_S << " Input and " << N_Group_T << " AdEx Neurons for " << stime << " seconds\n";

    /// Initializing AdEx Neurons
    std::vector<double> init_vtrest(N_Group_T, vtrest);
    std::vector<double> init_vm(N_Group_T, vtrest+ 0.005);
    std::vector<double> init_I(N_Group_T, 0);
    std::vector<double> init_x(N_Group_T, 0);
    AdEx AdEx_neurons(init_vtrest, init_vm, init_I, init_x);

    /// Initializing Input Neurons firing randomly
    std::vector<double> times;
    std::vector<int> spike_ids;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib_neurons(0., 1.);
    const int timesteps = (int)(stime/defaultclock_dt);
    for(auto t = 0; t < timesteps; t++) {
        for (auto i = 0; i < N_Group_S; i++) {
            if (distrib_neurons(gen) < 0.3) {
                times.push_back(t*defaultclock_dt);
                spike_ids.push_back(i);
            }
        }
    }
    Inputs input_neurons(N_Group_S, spike_ids, times);

    /// Initializing synapses between Input and AdEx Neurons
    std::vector<double> init_input_FBp(N_Group_S*N_Group_T, 0);
    std::vector<double> init_input_FBn(N_Group_S*N_Group_T, 0);
    std::vector<double> init_input_R(N_Group_S*N_Group_T, 1);
    std::vector<double> init_input_U(N_Group_S*N_Group_T);
    std::vector<double> init_input_A(N_Group_S*N_Group_T);
    const double input1_pos = 25;
    const double rad = 5;
    int init_const = 0;
    int index = 0;
    std::uniform_int_distribution<> distrib_synapses(0, 500);
    for(int i = 0; i < N_Group_S; i++){
        for(int j = 0; j < N_Group_T; j++){
            index = i*N_Group_T + j;
            init_const = distrib_synapses(gen);
            init_input_U[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;
            init_input_A[index] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;
        }
    }
    Synapses Input_AdEx_synapses(input_neurons, AdEx_neurons, init_input_FBp, init_input_FBn, init_input_R, init_input_U, init_input_A);

    /// Main simulation loop
    std::cout << "timesteps=" << timesteps << "\n";
    for(int t = 0; t < timesteps; t++){
        input_neurons.generate_spikes(t*defaultclock_dt);
        AdEx_neurons.solve_neurons();
        input_neurons.print_spikes();
        AdEx_neurons.print_spikes();
        Input_AdEx_synapses.UpdateSynapses_pre(t*defaultclock_dt);
        Input_AdEx_synapses.UpdateSynapses_post(t*defaultclock_dt);
    }
    /// Normally around 3.6MHz for 100x1000, 0.01 (s)
    std::cout << "Spike frequency is: " << (double)AdEx_neurons.get_total_spikes() / stime << "Hz" << std::endl;

	return 0;
}
