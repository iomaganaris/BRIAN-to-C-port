/**
 * @file neuron.h
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2014
 * @brief Header file containing all the declarations of Neuron classes and neurons' functions.
 */


/**
 * @brief Declaration of struct Neuron, that contains all the variables
 * of a neuron.
 *
 */
typedef struct {
    double vt;	/**< Voltage threshold. */
    double vm;	/**< Membrane potential.	*/
    double I;	/**< Neuron input current.	*/
    double x;	/**< Adaption variable (w).	*/
    int Spike;	/**< Variable to show if the neuron has spiked in a specific moment.	*/
} Neuron;
/**
 * @brief Declaration of struct Poisson, that contains all the variables
 * of a dummy input Poisson neuron that produces spikes randomly based in 
 * a poisson distribution.
 *
 */
typedef struct {
	double GaussArray;	/**< Gauss initial value. */
    int Spike;	/**< Variable to show if the neuron has spiked in a specific moment. */
} Poisson;
/**
 * @brief Function to print all the variables of the neurons.
 *
 * @param neurons Array of dummy and real (ADEX) neurons.
 * @param N Number of neurons to print (dummy or ADEX).
 */
void print_neurons(Neuron* neurons, int N);
/**
 * @brief Function to reset all the variables of a neuron if it has spiked.
 *
 * @param neuron Neuron to reset.
 */
void resetNeuron(Neuron* neuron);
/**
 * @brief Function to calculate which Poisson neurons have spiked in a specific moment.
 *
 * @param input Array containing all the poisson neurons.
 * @param N_S Number of dummy-input (source) neurons.
 * @param N_T Number of ADEX (target) neurons.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 */
void PoissonThreshold(Poisson* input, int N_S, int N_T, int* SpikeArray);
/**
 * @brief Function to update the variables of the neurons. It is called in every timestep.
 *
 * @param neurons Array containing all the neurons.
 * @param N Number of neurons to update.
 * @param SpikeArray Array that contains if a neuron had spiked in a specific moment.
 */
void SolveNeurons(Neuron* neurons, int N, int *SpikeArray);
