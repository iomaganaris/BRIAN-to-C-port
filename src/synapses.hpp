/**
 * @file synapses.h
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2014
 * @brief Header file containing all the declarations of Synapse class and synapses' functions.
 */
#include <vector>

#include "neuron.hpp"

/**
 * @brief Declaration of struct Synapse, that contains all the variables
 * of a synapse.
 *
 */
class Synapses {
    std::vector<double> w;	/**< Weight of a synapse. (Was present in BRIAN but never used for the ADEX with STDP simulation.) */
    std::vector<double> FFp;	/**< FFp variable of a synapse. (x+) */
    std::vector<double> FBp;	/**< FBp variable of a synapse. (y+) */
    std::vector<double> FBn;	/**< FBn variable of a synapse. (y-) */
    std::vector<double> R;	/**< R value of a synapse. (r)	*/
    std::vector<double> u;	/**< u value of a synapse. (p)	*/
    std::vector<double> U;	/**< U value of a synapse. (P)	*/
    std::vector<double> A;	/**< A value of a synapse. (q)	*/
    std::vector<double> lastupdate;	/**< Last time a synapse was updated	*/
    std::vector<double> target_I;	/**< The I value for the postsynaptic neuron. */
    int N_Group_S;  /**< The number of pre synaptic neurons. */
    int N_Group_T;  /**< The number of post synaptic neurons. */
    Neurons pre_neurons;
    AdEx post_neurons;
public:
/**
 * Constructor
 */
Synapses(Neurons& pre_neurons, AdEx& post_neurons): pre_neurons(pre_neurons), post_neurons(post_neurons) {
        N_Group_S = pre_neurons.get_n_neurons();
        N_Group_T = post_neurons.get_n_neurons();
        init();
}

/**
 * @brief Initialize the synapses vectors
 */
void init();

/**
 * @brief Function to print all the variables of the synapses.
 */
void print_synapses();
/**
 * @brief Function that encapsulates the presynaptic expression of STDP
 * and updates all the variables of a synapse that had a presynaptic event.
 *
 * @param Synapses 2D matrix containing all the synapses.
 * @param neurons Array of dummy and real (ADEX) neurons.
 * @param N_S Number of dummy-input (source) neurons.
 * @param N_Group_S Number of ADEX neurons that are connected with other ADEX neurons.
 * @param N_Group_T Number of ADEX (target) neurons that are connected with other ADEX or input neurons.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 * @param t Time.
 */
void UpdateSynapses_pre(double t);
/**
 * @brief Function that encapsulates the postsynaptic expression of STDP
 * and updates all the variables of a synapse that had a postsynaptic event.
 *
 * @param Synapses 2D matrix containing all the synapses.
 * @param N_S Number of dummy-input (source) neurons.
 * @param N_Group_S Number of ADEX neurons that are connected with other ADEX neurons.
 * @param N_Group_T Number of ADEX (target) neurons that are connected with other ADEX or input neurons.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 * @param t Time.
 */
void UpdateSynapses_post(double t);
};


