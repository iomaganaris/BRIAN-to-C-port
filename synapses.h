#include "neuron.h"

typedef struct {
	int conn;
	double w;		// if w = 0 no connection maybe (?)
	double FFp;
	double FBp;
	double FBn;
	double R;
	double u;
	double U;
	double A;
	double lastupdate;
	double target_I;
} Synapse;

void print_synapses(Synapse** syn, int N_S, int N_T);
void UpdateSynapses_pre(Synapse** Synapses, Neuron* neurons, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t);
void UpdateSynapses_post(Synapse** Synapses, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t);
