/**
 * @file synapses.h
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2014
 * @brief Header file containing all the declarations of Synapse class and synapses' functions.
 */
#include "neuron.h"

/**
 * @brief Declaration of struct Synapse, that contains all the variables
 * of a synapse.
 *
 */
typedef struct {
	int conn;	/**< Variable that expresses if a given synapse between two neurons exists. */
	double w;	/**< Weight of a synapse. (Was present in BRIAN but never used for the ADEX with STDP simulation.) */
	double FFp;	/**< FFp variable of a synapse. (x+) */
	double FBp;	/**< FBp variable of a synapse. (y+) */
	double FBn;	/**< FBn variable of a synapse. (y-) */
	double R;	/**< R value of a synapse. (r)	*/
	double u;	/**< u value of a synapse. (p)	*/
	double U;	/**< U value of a synapse. (P)	*/
	double A;	/**< A value of a synapse. (q)	*/
	double lastupdate;	/**< Last time a synapse was updated	*/
	double target_I;	/**< The I value for the postsynaptic neuron. */
} Synapse;

/**
 * @brief Function to print all the variables of the synapses.
 *
 * @param syn 2D matrix containing all the synapses.
 * @param N_S Number of dummy-input (source) neurons.
 * @param N_T Number of ADEX (target) neurons.
 */
void print_synapses(Synapse** syn, int N_S, int N_T);
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
void UpdateSynapses_pre(Synapse** Synapses, Neuron* neurons, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t);
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
void UpdateSynapses_post(Synapse** Synapses, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t);
