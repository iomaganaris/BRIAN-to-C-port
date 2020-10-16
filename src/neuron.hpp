/**
 * @file neuron.h
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2014
 * @brief Header file containing all the declarations of Neuron classes and neurons' functions.
 */

#include <vector>

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
    spikes.resize(n_neurons);
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
/**
 * @brief Function to update the variables of the neurons. It is called in every timestep.
 *
 * @param neurons Array containing all the neurons.
 * @param N Number of neurons to update.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 */
virtual void SolveNeurons() = 0;
};

class Input : public Neurons {
    std::vector<int> spike_replay;
public:
Input(const int nNeurons, std::vector<int> &ids, std::vector<int> &times) : Neurons(ids.size()) {
    spike_replay.resize(ids.size()*times.size());
}
void SolveNeurons() override;
};

class AdEx : public Neurons {
    std::vector<double> vt;	/**< Voltage threshold. */
    std::vector<double> vm;	/**< Membrane potential.	*/
    std::vector<double> I;	/**< Neuron input current.	*/
    std::vector<double> x;	/**< Adaption variable (w).	*/
public:
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
void SolveNeurons() override;
};

/**
 * @brief Declaration of struct Poisson, that contains all the variables
 * of a dummy input Poisson neuron that produces spikes randomly based in
 * a poisson distribution.
 *
 */
class PoissonNeurons: public Neurons {
    std::vector<double> GaussArray;	/**< Gauss initial value. */
public:
/**
 * @brief Function to calculate which Poisson neurons have spiked in a specific moment.
 *
 * @param input Array containing all the poisson neurons.
 * @param N_S Number of dummy-input (source) neurons.
 * @param N_T Number of ADEX (target) neurons.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 */
    void SolveNeurons() override;
};
