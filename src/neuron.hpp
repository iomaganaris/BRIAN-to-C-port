/**
 * @file neuron.hpp
 * @author Ioannis Magkanaris
 * @date 4 November 2020
 * @brief Header file containing all the declarations of Neuron classes and neurons' functions
 */

#include <vector>

/**
 * @brief Declaration of parent Neuron class.
 *
 */
class Neurons {
protected:
    std::vector<int> spikes; /// Vector of spikes
    int n_neurons;           /// Number of neurons
    int n_spikes;            /// Total number of spikes
public:

/**
 * @brief Constructor for Neurons
 *
 * @param n_neurons Number of neurons
 */
explicit Neurons(const int n_neurons): n_neurons(n_neurons){
    spikes.resize(n_neurons, 0);
}

/**
 * @brief Get number of neurons
 *
 */
int get_n_neurons() const {
    return n_neurons;
}

/**
 * @brief Get total spikes of the simulation
 *
 */
int get_total_spikes() const {
    return n_spikes;
}

/**
 * @brief Return the spikes at the current timestep
 *
 */
const std::vector<int>& get_spikes() {
    return spikes;
}
/**
 * @brief Function to print all the variables of the neurons.
 *
 */
void print_spikes();
};


/**
 * @brief Declaration of Input Neurons class.
 *
 * The Input Neurons have a vector for the id of the Input Neuron which generates a spike and a
 * corresponging vector of the timings that these spikes are generated. The two vectors should be
 * initialized sorted by the spike_times.
 *
 */
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
Inputs(const int n_input_neurons, std::vector<int> &ids, std::vector<double> &times) :
    Neurons(n_input_neurons), spike_ids(ids), spike_times(times) {};

/**
 * @brief Generates the spikes based on spike_ides and spike_times
 *
 * @param t Time of the simulation
 */
void generate_spikes(double t);

};


/**
 * @brief Declaration of AdEx Neurons class.
 * 
 * The behavior of AdEx neurons is described by the following differential equations:
 * eqs_neuron = """
 * dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-vt)/DeltaT)+I-x)/C : volt
 * dvt/dt=-(vt-vtrest)/tauvt : volt
 * dx/dt=(c*(vm-EL)-x)/tauw : amp #In the standard formulation x is w
 * I : amp
 *
 */
class AdEx : public Neurons {
    std::vector<double> vt;	/// Voltage threshold
    std::vector<double> vm;	/// Membrane potential
    std::vector<double> I;	/// Neuron input current
    std::vector<double> x;	/// Adaption variable (w)
public:

/**
 * @brief Constructor for AdEx neurons
 *
 * @param init_vt The vt initial vector
 * @param init_vm The vm initial vector
 * @param init_I The I initial vector
 * @param init_x The x initial vector
 * 
 */
AdEx(const std::vector<double>& init_vt,
    const std::vector<double>& init_vm,
    const std::vector<double>& init_I,
    const std::vector<double>& init_x) :
    Neurons(init_vt.size()),
        vt(init_vt),
        vm(init_vm),
        I(init_I),
        x(init_x) {};

/**
 * @brief Function to print all the variables of the neurons.
 *
 * @param pre_spikes Vector of spikes of pre-synaptic neurons
 * @param synapses_I Vector of current for all the synapses between pre- and post-synaptic neurons
 */
void update_I(std::vector<int>& pre_spikes, std::vector<double>& synapses_I);

/**
 * @brief Function to print all the variables of the neurons.
 *
 */
void print_neurons() const;

/**
 * @brief Function to reset all the variables of a neuron if it has spiked.
 *
 * @param id Neuron to reset.
 */
inline void resetNeuron(int id);

/**
 * @brief Function that solves the equations of the neurons and updates the spike vector.
 *
 */
void solve_neurons();

};
