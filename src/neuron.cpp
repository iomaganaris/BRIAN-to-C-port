/**
 * @file neuron.cpp
 * @author Ioannis Magkanaris
 * @date 4 November 2020
 * @brief File containing all the functions of Neurons
 *
 * It contains the declarations of all the types of Neurons
 */

#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "neuron.hpp"

void Neurons::print_spikes() {
    std::cout << "\nSpikes\n";
    for(const auto& spikes_local : spikes){
        std::cout << spikes_local << ", ";
    }
    std::cout << std::endl;
}

void Inputs::generate_spikes(const double t) {
    std::fill(spikes.begin(), spikes.end(), 0);
    while(spike_times[current_index] == t) {
        spikes[spike_ids[current_index]] = 1;
        current_index++;
    }
}

void AdEx::print_neurons() const{
	std::cout << "\nvt\n";
	for(const auto& vt_local : vt){
		std::cout << vt_local << ", ";
	}
    std::cout << "\nvm\n";
    for(const auto& vm_local : vm){
        std::cout << vm_local << ", ";
    }
    std::cout << "\nI\n";
    for(const auto& I_local : I){
        std::cout << I_local << ", ";
    }
    std::cout << "\nx\n";
    for(const auto& x_local : x){
        std::cout << x_local << ", ";
    }
    std::cout << "\nSpikes\n";
    for(const auto& spikes_local : spikes){
        std::cout << spikes_local << ", ";
    }
	std::cout << std::endl;
}

inline void AdEx::resetNeuron(const int id) {
    vm[id] = Vr;
    x[id] += b;
    vt[id] = VTmax;
}

void AdEx::solve_neurons() {
    for(auto id = 0; id < n_neurons; id++){
        double _vm, _vt, _x;
        _vm = (gL*(EL-vm[id])+gL*DeltaT*exp((vm[id]-vt[id])/DeltaT)+I[id]-x[id])/C;
        _vt = -(vt[id]-vtrest)/tauvt;
        _x = (c*(vm[id]-EL)-x[id])/tauw;
        vm[id] += _vm * defaultclock_dt;
        vt[id] += _vt * defaultclock_dt;
        x[id] += _x * defaultclock_dt;
        if(vm[id] > vt[id]){
            resetNeuron(id);
            spikes[id] = 1;
            n_spikes++;
        }
        else spikes[id] = 0;
    }
}

void AdEx::update_I(std::vector<int>& pre_spikes, std::vector<double>& synapses_I) {
    const auto n_pre_neurons = pre_spikes.size();
    const auto n_post_neurons = n_neurons;
    auto index = 0;
    for(auto i = 0; i < n_pre_neurons; i++) {
        for(auto j = 0; j < n_post_neurons; j++) {
            if(pre_spikes[i]) {
                index = i * n_post_neurons + j;
                I[j] = synapses_I[index];
            }
        }
    }
}
