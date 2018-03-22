typedef struct {
    double vt;// = vtrest;
    double vm;// = EL;
    double I;// = 0;
    double x;// = 0;
    int Spike;// = false;
} Neuron;

typedef struct {
	double GaussArray;
    int Spike; 
} Poisson;

void print_neurons(Neuron* neurons, int N);
void resetNeuron(Neuron* neuron);
void PoissonThreshold(Poisson* input, int N_S, int N_T, int* SpikeArray);
void SolveNeurons(Neuron* neurons, int N, int *SpikeArray);
