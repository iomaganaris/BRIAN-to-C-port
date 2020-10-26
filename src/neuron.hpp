/**
 * @file neuron.h
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2014
 * @brief Header file containing all the declarations of Neuron classes and neurons' functions.
 */

#include <vector>

#include "util.hpp"

/**
 * @brief Declaration of parent Neuron class.
 *
 */
class Neurons {
protected:
    std::vector<int> spikes;	/**< Variable to show if the neuron has spiked in a specific moment.	*/
    int n_neurons;
public:
int get_n_neurons() const {
    return n_neurons;
}
/**
 * Constructor
 */
explicit Neurons(const int n_neurons): n_neurons(n_neurons){
    spikes.resize(n_neurons, 0);
}
const std::vector<int>& get_spikes() {
    return spikes;
}
/**
 * @brief Function to print all the variables of the neurons.
 *
 * @param neurons Array of dummy and real (ADEX) neurons.
 * @param N Number of neurons to print (dummy or ADEX).
 */
void print_spikes();
};

class Inputs : public Neurons {
    std::vector<int> spike_ids;
    std::vector<double> spike_times;
    int current_index = 0;
public:
/**
 * @brief Constructor for input spikes
 *
 * @param n_input_neurons Number of input neurons
 * @param ids IDs of the neurons that spike
 * @param times The corresponding timing for every ID
 */
Inputs(const int n_input_neurons, std::vector<int> &ids, std::vector<double> &times) : Neurons(n_input_neurons), spike_ids(ids), spike_times(times) {
    /// Sort spikes by times to make it easier to read later
    //sort_second(spike_ids, spike_times);
};
void generate_spikes(double t);
};

class AdEx : public Neurons {
    std::vector<double> vt;	/**< Voltage threshold. */
    std::vector<double> vm;	/**< Membrane potential.	*/
    std::vector<double> I;	/**< Neuron input current.	*/
    std::vector<double> x;	/**< Adaption variable (w).	*/
public:
AdEx(const std::vector<double>& init_vt, const std::vector<double>& init_vm, const std::vector<double>& init_I, const std::vector<double>& init_x) : Neurons(init_vt.size()), vt(init_vt), vm(init_vm), I(init_I), x(init_x) {};
void update_I(std::vector<int>& pre_spikes, std::vector<double>& synapses_I);
/**
 * @brief Function to print all the variables of the neurons.
 *
 * @param neurons Array of dummy and real (ADEX) neurons.
 * @param N Number of neurons to print (dummy or ADEX).
 */
void print_neurons() const;
/**
 * @brief Function to reset all the variables of a neuron if it has spiked.
 *
 * @param neuron Neuron to reset.
 */
inline void resetNeuron(int id);
void solve_neurons();
};

/**
 * @brief Declaration of struct Poisson, that contains all the variables
 * of a dummy input Poisson neuron that produces spikes randomly based in
 * a poisson distribution.
 *
 */
//class PoissonNeurons: public Neurons {
//    std::vector<double> GaussArray;	/**< Gauss initial value. */
//public:
/**
 * @brief Function to calculate which Poisson neurons have spiked in a specific moment.
 *
 * @param input Array containing all the poisson neurons.
 * @param N_S Number of dummy-input (source) neurons.
 * @param N_T Number of ADEX (target) neurons.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 */
//    void SolveNeurons() override;
//};
