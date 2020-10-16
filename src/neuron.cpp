#include <cmath>
#include <cstdlib>
#include <iostream>

#include "constants.hpp"
#include "neuron.hpp"

void Neurons::print_spikes() {
    std::cout << "\nSpikes\n";
    for(const auto& spikes_local : spikes){
        std::cout << spikes_local < ", ";
    }
    std::cout << std::endl;
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
        std::cout << spikes_local < ", ";
    }
	std::cout << std::endl;
}

inline void AdEx::resetNeuron(const int id) {
    vm[id] = Vr;
    x[id] += b;
    vt[id] = VTmax;
}

void AdEx::SolveNeurons(){
    for(auto id = 0; id < n_neurons; id++){
        double _vm, _vt, _x;
        _vm = (gL*(EL-vm[id])+gL*DeltaT*exp((vm[id]-vt[id])/DeltaT)+I[id]-x[id])/C;
        //printf("_vm = %.20f\n", _vm);
        _vt = -(vt[id]-vtrest)/tauvt;
        //printf("_vt = %.20f\n", _vm);
        _x = (c*(vm[id]-EL)-x[id])/tauw;
        //printf("_x = %.20f\n", _x);
        vm[id] += _vm * defaultclock_dt;
        vt[id] += _vt * defaultclock_dt;
        x[id] += _x * defaultclock_dt;
        if(vm[id] > vt[id]){
            //printf("Reset\n");
            resetNeuron(id);
            spikes[id] = 1;
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

/*double random_0_1(){
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
}*/
