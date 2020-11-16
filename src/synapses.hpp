/**
 * @file synapses.hpp
 * @author Ioannis Magkanaris
 * @date 4 November 2020
 * @brief Header file containing all the declarations of Synapse class and synapses' functions
 */

#include <vector>

#include "neuron.hpp"

/**
 * @brief Declaration of struct Synapse, that contains all the variables
 * of a synapse.
 *
 */
class Synapses {
    std::vector<double> w;	         /// Weight of a synapse
    std::vector<double> FFp;	     /// FFp variable of a synapse (x+)
    std::vector<double> FBp;	     /// FBp variable of a synapse (y+)
    std::vector<double> FBn;	     /// FBn variable of a synapse (y-)
    std::vector<double> R;	         /// R value of a synapse (r)
    std::vector<double> u;	         /// u value of a synapse (p)
    std::vector<double> U;	         /// U value of a synapse (P)
    std::vector<double> A;	         /// A value of a synapse (q)
    std::vector<double> lastupdate;	 /// Last time a synapse was updated
    std::vector<double> target_I;	 /// The I value for the postsynaptic neuron
    int N_Group_S;                   /// The number of pre synaptic neurons
    int N_Group_T;                   /// The number of post synaptic neurons
    Neurons& pre_neurons;            /// The pre-synaptic neurons
    AdEx& post_neurons;              /// The post-synaptic neurons
public:
/**
 * @brief Constructor of Synapses
 *
 * Creates the connectivity map of the Synapses as SoA
 *
 * @param pre_neurons The pre-synaptic Neurons
 * @param post_neurons The post-synaptic Neurons
 * @param init_FBp The initial FBp vector
 * @param init_FBn The initial FBn vector
 * @param init_R The initial R vector
 * @param init_U The initial U vector
 * @param init_A The initial A vector
 *
 */
Synapses(Neurons& pre_neurons,
    AdEx& post_neurons,
    const std::vector<double>& init_FBp,
    const std::vector<double>& init_FBn,
    const std::vector<double>& init_R,
    const std::vector<double>& init_U,
    const std::vector<double>& init_A) :
        pre_neurons(pre_neurons),
        post_neurons(post_neurons),
        FBp(init_FBp),
        FBn(init_FBn),
        R(init_R),
        U(init_U),
        A(init_A) {
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
 *
 */
void print_synapses();

/**
 * @brief Function that encapsulates the presynaptic expression of STDP
 * and updates all the variables of a synapse that had a presynaptic event.
 *
 * @param t Simulation time
 *
 */
void UpdateSynapses_pre(double t);
/**
 * @brief Function that encapsulates the postsynaptic expression of STDP
 * and updates all the variables of a synapse that had a postsynaptic event.
 *
 * @param t Simulation time
 *
 */
void UpdateSynapses_post(double t);
};


