#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//# Adex Parameters
double C = 281;	//pF
double gL = 30; //nS
double taum = 281 / 30;	// C/gL	// double initilization of taum(?)
double EL = -70.6;	//mV
double DeltaT = 2;	//mV
double vti = -50.4;	//mV
//#vtrest = vti + 5 * DeltaT
double vtrest = -45;	//mV
double VTmax = 18;	//mV
double tauvt = 50;	//ms

double tauw = 144;	//ms
double c = 4;	//ns
double b = 0.805;	//nA
double Vr = -70.6;	//mV

typedef struct {
    double vt;// = vtrest;
    double vm;// = EL;
    double I;// = 0;
    double x;// = 0;
    int Spike;// = false;
} Neuron;
void resetNeuron(Neuron* neuron) {
    neuron->vm = Vr;
    neuron->x += b;
    neuron->vt = VTmax;
}
typedef struct {
	double* GaussArray;
    int Spike; 
} Poisson;
void SolveNeurons(Neuron* neurons, int N, int *SpikeArray){
	for(int i = 0; i < N; i++){
		if (rand() % 10 < 5) {
			SpikeArray[i] = 0;
		}	
		else {
			SpikeArray[i] = 1;
		}
	}

}

int main(void){

	int nruns = 1;

	for(int nrun = 0; nrun < nruns; nrun++){
		double realtime = 0;
		double stime = 1; //second
		double stime2 = 50; //second

		double resolution_export = 10; //every x ms

		int N = 10;
		double taum = 10; //ms
		double Ee = 0; //mV
		double tuae = 2; //ms
		double Fon = 50; //Hz
		double Foff = 3; //Hz

		double s = 100*1e-10; 
		double Amax = 2.0;
		double Amin = 0;
		double Ainit = 0.1;
		double Umax = 1.0;
		double Umin = 0;
		double Uinit = 0.1;

		double dFBn = 0;
		double dFBp = 0;
		double dFFp = 0;

		//#Short-term plasticity params
	    double tau_u = 50;	//ms
	    double tau_r = 200;	//ms

	    //#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
	    double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
	    double AFBn = params[0];
	    double tau_FBn = params[1]*1e3;	//ms
	    double AFBp = params[2];
	    double tau_FBp = params[3]*1e3;	//ms
	    double AFFp = params[4];
	    double tau_FFp = params[5]*1e3;	//ms
	    //#etaU = 0.35
	    double etaU = 0.15;
	    double etaA = 0.15;
	    //#etaA = 0.35    

        double defaultclock_dt = 1;	//ms

	    
	    /*//# Adex Parameters
	    double C = 281;	//pF
	    double gL = 30; //nS
	    double taum = C / gL;		// double initilization of taum(?)
	    double EL = -70.6;	//mV
	    double DeltaT = 2;	//mV
	    double vti = -50.4;	//mV
	    //#vtrest = vti + 5 * DeltaT
	    double vtrest = -45;	//mV
	    double VTmax = 18;	//mV
	    double tauvt = 50;	//ms

	    double tauw = 144;	//ms
	    double c = 4;	//ns
	    double b = 0.805;	//nA
	    double Vr = -70.6;	//mV*/
	    
	    /*eqs_neuron = """
        dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-vt)/DeltaT)+I-x)/C : volt
        dvt/dt=-(vt-vtrest)/tauvt : volt
        dx/dt=(c*(vm-EL)-x)/tauw : amp #In the standard formulation x is w
        I : amp
    	"""*/

    	double input1_pos = 25;
    	double input2_pos = 75;
		double rad = 5;

	    //Define input 1
	    double *F_input1;
	    F_input1 = (double*)malloc(sizeof(double)*N);
	    //F_input1[input1_pos-rad:input1_pos+rad] = Fon
	    for(int i = 0; i < N; i++){
	    	F_input1[i] = Foff;			//maybe is not needed
	    	F_input1[i] = exp(-(pow((i+1)-input1_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    }

	    //Define input 2
	    double *F_input2;
	    F_input2 = (double*)malloc(sizeof(double)*N);
	    //F_input2[input2_pos-rad:input2_pos+rad] = Fon
	    for(int i = 0; i < N; i++){
	    	F_input2[i] = Foff;			//maybe is not needed
	    	F_input2[i] = exp(-(pow((i+1)-input2_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    }

	    int N_S = N;
	    Poisson *input;
	    input = (Poisson*)malloc(sizeof(Poisson)*N_S);
	    for(int i = 0; i<N_S; i++){				// Initialization of Poisson Neurons
	    	input[i].GaussArray = F_input1;
	    	input[i].Spike = 0;
	    }
	    int N_T = N;	//for the simulation we have, normaly is N
	    Neuron *neurons;
	    neurons = (Neuron*)malloc(sizeof(Neuron)*N_T);
	    for(int i = 0; i<N_T; i++){				// Initiliazation of Neurons
	    	neurons[i].vt = vtrest;
	    	neurons[i].vm = EL;
	    	neurons[i].I = 0;
	    	neurons[i].x = 0;
	    	neurons[i].Spike = 0;
	    }
	    printf("neurons->vm\n");
	    for(int i=0; i<N_T; i++){
	    	printf("%lf\n",neurons[i].vm);
	    }

	    int *SpikeArray;
	    SpikeArray = (int*)malloc(sizeof(int)*N);

	    //To be initialized
	    /*neurons.vt = vtrest
    	neurons.vm = EL
    	neurons.I = 0
    	neurons.x = 0*/

    	/*model='''w : 1
             FFp : 1
             FBp : 1
             FBn : 1
             R : 1
             u : 1
             U : 1
             A : 1         
             dFFp/dt=-FFp/tau_FFp : 1 (event-driven)
             dFBp/dt=-FBp/tau_FBp : 1 (event-driven)
             dFBn/dt=-FBn/tau_FBn : 1 (event-driven)
             dR/dt=(1-R)/tau_r : 1 (event-driven)
             du/dt=(U-u)/tau_u : 1 (event-driven)            
             '''*/
	    //Synapse syn[N_S][N_T] = CreateSynapses(input, neurons);	//2D Array of synapses. Each element/synapse has it's variables embedded.
	    												//Must find a way to describe how they are connected
	    //syn.connect_one_to_one(input, neurons)
	    //syn[:,:]=True	//what is that?	#TO DO
	    //syn.FBp=0		#TO DO
	    //syn.FBn=0		#TO DO
	    //syn.R=1		#TO DO
	    //syn.U='rand()*Uinit'
	   	//syn.A='rand()*Ainit'
	    //syn.U[:]=Umin
	    //syn.U[:]=0.5
	    
	    /*for(int i =0; i < N; i++){					//Define gaussian input
	    	for(int j = 0; j < N; j++){ 
	        syn[i][j].U = exp(-(((pow((i+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;
	        syn[i][j].A = exp(-(((pow((i+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;
	    	}
		}*/

	    srand(time(NULL));

		int timesteps = stime*100/resolution_export;
		for(int t = 0; t < timesteps; t++){			//add monitors for the variables we care about
			SolveNeurons(neurons, N_T, SpikeArray);
			printf("Loop %d\n",t);
			for(int i = 0; i < N_T; i++){
				printf("%d\n",SpikeArray[i]);
			}
			//UpdateSynapses_pre(SpikeArray);
			//UpdateSynapses_post(SpikeArray);
		}
	}
	return 0;
}
