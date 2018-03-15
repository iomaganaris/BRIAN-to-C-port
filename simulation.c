#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double defaultclock_dt = 1*1e-3;	//ms

//double taum_S = 10 * 1e-3; //ms
double Ee = 0 * 1e-3; //mV
double tuae = 2 * 1e-3; //ms
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
double tau_u = 50 * 1e-3;	//ms
double tau_r = 200 * 1e-3;	//ms

//#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
//double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
double AFBn = 0.1771;
double tau_FBn = 0.0327 * 1e3 * 1e-3;	//ms
double AFBp = 0.1548;
double tau_FBp = 0.2302 * 1e3 * 1e-3;	//ms
double AFFp = 0.0618;
double tau_FFp = 0.0666 * 1e3 * 1e-3;	//ms
//#etaU = 0.35
double etaU = 0.15;
double etaA = 0.15;
//#etaA = 0.35    


//# Adex Parameters
double C = 281*1e-12;	//pF
double gL = 30*1e-9; //nS
double taum = 281*1e-12 / 30*1e-9;	// C/gL	// double initilization of taum(?)
double EL = -70.6*1e-3;	//mV
double DeltaT = 2*1e-3;	//mV
double vti = -50.4*1e-3;	//mV
//#vtrest = vti + 5 * DeltaT
double vtrest = -45*1e-3;	//mV
double VTmax = 18*1e-3;	//mV
double tauvt = 50*1e-3;	//ms

double tauw = 144*1e-3;	//ms
double c = 4*1e-9;	//ns
double b = 0.0805*1e-9;	//nA
double Vr = -70.6*1e-3;	//mV

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

void resetNeuron(Neuron* neuron) {
    neuron->vm = Vr;
    neuron->x += b;
    neuron->vt = VTmax;
}

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
double random_0_1(){
	return (double)rand()/(double)((unsigned)RAND_MAX+1);
}
void PoissonThreshold(Poisson* input, int N_S, int N_T, int* SpikeArray){
	for(int i = 0; i<N_S; i++){
		if(random_0_1()<input[i].GaussArray*defaultclock_dt){
			input[i].Spike = 1;
			SpikeArray[N_T+i] = 1;
			printf("PoissonThreshold i= %d\n",i);
		} 
		else{
			input[i].Spike = 0;
			SpikeArray[N_T+i] = 0;
		}
	}
	return;
}

void SolveNeurons(Neuron* neurons, int N, int *SpikeArray){
    for(int i = 0; i < N; i++){
        /*if (rand() % 10 < 5) {
            SpikeArray[i] = 0;
        }
        else {
            SpikeArray[i] = 1;
        }*/
        double _vm, _vt, _x;
        _vm = (gL*(EL-neurons[i].vm)+gL*DeltaT*exp((neurons[i].vm-neurons[i].vt)/DeltaT)+neurons[i].I-neurons[i].x)/C;
        //printf("_vm = %.20f\n", _vm);
        _vt = -(neurons[i].vt-vtrest)/tauvt;
        //printf("_vt = %.20f\n", _vm);
        _x = (c*(neurons[i].vm-EL)-neurons[i].x)/tauw;
        //printf("_x = %.20f\n", _x);
        neurons[i].vm += _vm * defaultclock_dt;
        neurons[i].vt += _vt * defaultclock_dt;
        neurons[i].x += _x * defaultclock_dt;
        if(neurons[i].vm > neurons[i].vt){
            //printf("Reset\n");
            resetNeuron(&neurons[i]);
            SpikeArray[i] = 1;
        }
        else SpikeArray[i] = 0;
    }
}

void UpdateSynapses_pre(Synapse** Synapses, Neuron* neurons, int N_S, int N_T, int* SpikeArray, double t){
	//an baloume ekswteriki eisodo thelei i<N+1 logika
	for (int i = 0; i < N_S; i++){
		if (SpikeArray[i] > 0){
			printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
			for (int j = 0; j < N_T; j++){
				//auto borei na ginei pio apodotika alla diskoleuei, gi arxi trwme xwro me ta struct synapses
				if (Synapses[i][j].conn){
					printf("i = %d, j = %d\n",i,j);
					Synapses[i][j].FFp = Synapses[i][j].FFp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FFp);
					Synapses[i][j].FBn = Synapses[i][j].FBn * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBn);
					Synapses[i][j].u = Synapses[i][j].U + (-Synapses[i][j].U + Synapses[i][j].u) * exp(-(-Synapses[i][j].lastupdate + t)/tau_u);
					Synapses[i][j].FBp = Synapses[i][j].FBp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBp);
					Synapses[i][j].R = (Synapses[i][j].R - 1) * exp(-(-Synapses[i][j].lastupdate + t)/tau_r) + 1;
					Synapses[i][j].target_I = s * Synapses[i][j].A * Synapses[i][j].R * Synapses[i][j].u;
					Synapses[i][j].U = Synapses[i][j].U + etaU * (-AFBn * Synapses[i][j].FBn * Synapses[i][j].FBp + AFBp * Synapses[i][j].FBp * Synapses[i][j].FFp);
					if (Synapses[i][j].U < Umin) Synapses[i][j].U = Umin;
					else if (Synapses[i][j].U > Umax) Synapses[i][j].U = Umax;
					Synapses[i][j].w = Synapses[i][j].U * Synapses[i][j].A;
					Synapses[i][j].FFp += 1;
					Synapses[i][j].R -= Synapses[i][j].R * Synapses[i][j].u;
					Synapses[i][j].u += Synapses[i][j].U * (1 - Synapses[i][j].u);
					Synapses[i][j].lastupdate = t;
				}
			}
		}
	}
	for (int j = 0; j < N_T; j++){
		for (int i = 0; i < N_S; i++){
			if (Synapses[i][j].conn && SpikeArray[i]){
				//printf("i: %d, j: %d\n",i,j);
				neurons[j].I = Synapses[i][j].target_I;
				printf("j: %d, neurons[j]= %lf",j,neurons[j].I);
				//break;
			}
		}
	}
}

void UpdateSynapses_post(Synapse** Synapses, int N_S, int N_T, int* SpikeArray, double t){
	
	for (int i = 0; i < N_T; i++){
	    if (SpikeArray[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
	            if (Synapses[j][i].conn){
	            	printf("i = %d, j = %d\n",i,j);
	            	fflush(stdout);
	                Synapses[j][i].FFp = Synapses[j][i].FFp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FFp);
	                Synapses[j][i].FBn = Synapses[j][i].FBn * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBn);
	                Synapses[j][i].u = Synapses[j][i].U + (-Synapses[j][i].U + Synapses[j][i].u) * exp(-(-Synapses[j][i].lastupdate + t)/tau_u);
	                Synapses[j][i].FBp = Synapses[j][i].FBp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBp);
	                Synapses[j][i].R = (Synapses[j][i].R - 1) * exp(-(-Synapses[j][i].lastupdate + t)/tau_r) + 1;
	                /*Synapses[j][i].A = Synapses[j][i].A + etaA * (AFFp * Synapses[j][i].FFp * Synapses[j][i].FBn);
	                double mean = 0;
	                int num = 0;
	                for (int k=0; k<N_S; k++){
	                    for (int l=0; l<N_T; l++){
	                        //auto edw to if boroume na to apofigoume an arxikopoioume tis metablites tou struct sto 0
	                        //if (Synapses[k][l].conn && (SpikeArray[k] || SpikeArray[l])){
	                    	if (Synapses[k][l].conn){
	                            mean += AFFp * Synapses[k][l].FFp * Synapses[k][l].FBn;
	                            num++;
	                        }
	                    }
	                }
	                mean = (double)mean / num;
	                printf("mean = %lf",mean);
	                Synapses[j][i].A = Synapses[j][i].A - etaA * 0.5 * mean; //amfibola swsto, sigoura mi apodotiko
	                if (Synapses[j][i].A < Amin) Synapses[j][i].A = Amin;
	                else if (Synapses[j][i].A > Amax) Synapses[j][i].A = Amax;
	                Synapses[j][i].w = Synapses[j][i].U * Synapses[j][i].A;
	                Synapses[j][i].FBp += 1;
	                Synapses[j][i].FBn += 1;*/

	            }
	        }
	    }
	}
	for (int i = 0; i < N_T; i++){
	    if (SpikeArray[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
	            if (Synapses[j][i].conn){
	            	Synapses[j][i].A = Synapses[j][i].A + etaA * (AFFp * Synapses[j][i].FFp * Synapses[j][i].FBn);
	                double mean = 0;
	                int num = 0;
	                for (int k=0; k<N_S; k++){
	                    for (int l=0; l<N_T; l++){
	                        //auto edw to if boroume na to apofigoume an arxikopoioume tis metablites tou struct sto 0
	                        //if (Synapses[k][l].conn && (SpikeArray[k] || SpikeArray[l])){
	                    	if (Synapses[k][l].conn){
	                            mean += AFFp * Synapses[k][l].FFp * Synapses[k][l].FBn;
	                            num++;
	                        }
	                    }
	                }
	                mean = (double)mean / num;
	                printf("mean = %lf",mean);
	                Synapses[j][i].A = Synapses[j][i].A - etaA * 0.5 * mean; //amfibola swsto, sigoura mi apodotiko
	                if (Synapses[j][i].A < Amin) Synapses[j][i].A = Amin;
	                else if (Synapses[j][i].A > Amax) Synapses[j][i].A = Amax;
	            }
	        }
	    }
	}
	for (int i = 0; i < N_T; i++){
	    if (SpikeArray[i] > 0){		// problima me ayto an exoyme allo group me neyrwnes san source kai target kai ton kommeno pinaka
	    	//printf("i= %d SpikeArray[i]= %d\n",i,SpikeArray[i]);
	    	//fflush(stdout);
	        for (int j = 0; j < N_S; j++){
	        	//printf("i = %d, j = %d Synapses[j][i].conn= %d\n",i,j,Synapses[j][i].conn);
				//fflush(stdout);
	            if (Synapses[j][i].conn){
	            	Synapses[j][i].w = Synapses[j][i].U * Synapses[j][i].A;
	                Synapses[j][i].FBp += 1;
	                Synapses[j][i].FBn += 1;
	                Synapses[j][i].lastupdate = t;
	            }
	        }
	    }
	}


}

void print_neurons(Neuron* neurons, int N){
	printf("vt\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].vt);
	}
	printf("\nvm\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].vm);
	}
	printf("\nI\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].I);
	}
	printf("\nx\n");
	for(int i = 0; i < N; i++){
		printf("%.8e, ",neurons[i].x);
	}
	printf("\nSpike\n");
	for(int i = 0; i < N; i++){
		printf("%d, ",neurons[i].Spike);
	}
	printf("\n");
}

void print_synapses(Synapse** syn, int N_S, int N_T){
	printf("conn\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			//printf("conn= %d, w= %lf, FFp= %lf, FBp= %lf, FBn= %lf, R= %lf, u= %lf, U= %lf, A= %lf, lastup= %lf, target_I= %lf\n",syn[i][j].conn,syn[i][j].w,syn[i][j].FFp,syn[i][j].FBp,syn[i][j].FBn,syn[i][j].R,syn[i][j].u,syn[i][j].U,syn[i][j].A,syn[i][j].lastupdate,syn[i][j].target_I);
			printf("%d ", syn[i][j].conn);
			if((i*N_T+j)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nw\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].w);
			if((i*N_T+j)%4 == 0) printf("\n");
		}
		//printf("\n");		
	}
	printf("\nFFp\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FFp);
		printf("\n");	
	}
	printf("\nFBp\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FBp);	
		printf("\n");
	}
	printf("\nFBn\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FBn);
		printf("\n");		
	}
	printf("\nR\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].R);
		printf("\n");
	}
	printf("\nu\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].u);
			if((i*N_T+j)%4 == 0) printf("\n");
		}
		//printf("\n");		
	}
	printf("\nU\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].U);
			if((i*N_T+j)%4 == 0) printf("\n");
		}
		//printf("\n");
	}
	printf("\nA\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e ", syn[i][j].A);	
			if((i*N_T+j)%4 == 0) printf("\n");
		}
		//printf("\n");
	}
	printf("\nlastupdate\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%lf, ", syn[i][j].lastupdate);		
		//printf("\n");
	}
	printf("\ntarget_I\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].target_I);
			if((i*N_T+j)%4 == 0) printf("\n");	
		}	
	}
	printf("\n");
}

int main(void){

	int nruns = 1;

	for(int nrun = 0; nrun < nruns; nrun++){
		double realtime = 0;
		double stime = 0.1; //second
		double stime2 = 50; //second

		double resolution_export = 10 * 1e-3; //every x ms

		int N = 1;
		int N_S = 100;
		int N_Group_S = N;
		int N_Group_T = N;	//for the simulation we have, normaly is N
		/*double taum = 10 * 1e-3; //ms
		double Ee = 0 * 1e-3; //mV
		double tuae = 2 * 1e-3; //ms
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
	    double tau_u = 50 * 1e-3;	//ms
	    double tau_r = 200 * 1e-3;	//ms

	    //#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
	    double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
	    double AFBn = params[0];
	    double tau_FBn = params[1]*1e3 * 1e-3;	//ms
	    double AFBp = params[2];
	    double tau_FBp = params[3]*1e3 * 1e-3;	//ms
	    double AFFp = params[4];
	    double tau_FFp = params[5]*1e3 * 1e-3;	//ms
	    //#etaU = 0.35
	    double etaU = 0.15;
	    double etaA = 0.15;
	    //#etaA = 0.35  */  

        //double defaultclock_dt = 1;	//ms

	    
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
	    F_input1 = (double*)malloc(sizeof(double)*N_S);
	    //F_input1[input1_pos-rad:input1_pos+rad] = Fon
	    printf("F_input1\n");
	    for(int i = 0; i < N_S; i++){
	    	F_input1[i] = Foff;			//maybe is not needed
	    	F_input1[i] = exp(-(pow((i+1)-input1_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    	printf("%lf ",F_input1[i]);
	    }
	    printf("\n");

	    //Define input 2
	    double *F_input2;
	    F_input2 = (double*)malloc(sizeof(double)*N_S);
	    //F_input2[input2_pos-rad:input2_pos+rad] = Fon
	    for(int i = 0; i < N_S; i++){
	    	F_input2[i] = Foff;			//maybe is not needed
	    	F_input2[i] = exp(-(pow((i+1)-input2_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    }

	    
	    Poisson *input;
	    input = (Poisson*)malloc(sizeof(Poisson)*N_S);
	    for(int i = 0; i<N_S; i++){				// Initialization of Poisson Neurons
	    	input[i].GaussArray = F_input1[i];
	    	input[i].Spike = 0;
	    }
	    
	    Neuron *neurons;
	    neurons = (Neuron*)malloc(sizeof(Neuron)*(N_Group_T+N_S));
	    for(int i = 0; i<N_Group_T; i++){				// Initiliazation of Neurons
	    	neurons[i].vt = vtrest;
	    	neurons[i].vm = EL; //vtrest + 0.005;//EL;
	    	neurons[i].I = 0;
	    	neurons[i].x = 0;
	    	neurons[i].Spike = 0;
	    }
	    /*printf("neurons->vm\n");
	    for(int i=0; i<N_T; i++){
	    	printf("%lf\n",neurons[i].vm);
	    }*/

	    int *SpikeArray;
	    SpikeArray = (int*)malloc(sizeof(int)*(N_Group_S+N_S));
	    
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

	    Synapse *syn[N_S+N_Group_S];
	    for(int i = 0; i < N_S+N_Group_S; i++) syn[i] = (Synapse*)malloc(sizeof(Synapse) * (N_Group_T+N_S));

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

	    // Initialization of Synapses for Neurons
	    for(int i = 0; i < N_Group_S; i++){
	    	for(int j = 0; j < N_Group_T; j++){
	    		syn[i][j].conn = 0;	// all connected
	    		syn[i][j].FBp = 0;
	    		syn[i][j].FBn = 0;
	    		syn[i][j].R = 1;
	    		syn[i][j].U = exp(-(((pow(((i*N_Group_T+j)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
	    		syn[i][j].A = exp(-(((pow(((i*N_Group_T+j)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
	    	}
	    }
	    // Initialization of Synapses for external input (bottom rows)
		for(int i = N_Group_S; i < N_Group_S+N_S; i++){
			for(int j = 0; j < N_Group_T + N_S; j++){
		    	syn[i][j].conn = 1;	// all connected
	    		syn[i][j].FBp = 0;
	    		syn[i][j].FBn = 0;
	    		syn[i][j].R = 1;
	    		syn[i][j].U = exp(-(((pow((((i-N_Group_S)*N_Group_T+j)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
	    		syn[i][j].A = exp(-(((pow((((i-N_Group_S)*N_Group_T+j)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
	    	}
		}
	    // Initialization of Synapses for external input ( right columns)
	    for(int i = 0; i < N_S + N_Group_S; i++){
		    for(int j = N_Group_T; j < N_Group_T+N_S; j++){
		    	syn[i][j].conn = 1;	// all connected
	    		syn[i][j].FBp = 0;
	    		syn[i][j].FBn = 0;
	    		syn[i][j].R = 1;
	    		syn[i][j].U = exp(-(((pow(((i*N_Group_T+j-N_Group_T)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
	    		syn[i][j].A = exp(-(((pow(((i*N_Group_T+j-N_Group_T)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
		    }
		}

	    /*printf("U:\n");
	    for(int i = 0 ; i < N_S+1; i++){
	    	for(int j = 0; j < N_T+1; j++){
	    		printf("%lf ",syn[i][j].U);
	    	}
	    	printf("\n");
	    }

	    printf("A:\n");
	    for(int i = 0 ; i < N_S+1; i++){
	    	for(int j = 0; j < N_T+1; j++){
	    		printf("%lf ",syn[i][j].A);
	    	}
	    	printf("\n");
	    }*/

	    srand(time(NULL));

		int timesteps = stime/defaultclock_dt;
		printf("timesteps=%d\n",timesteps);
		for(int t = 0; t < timesteps; t++){			//add monitors for the variables we care about
			printf("t: %.20lf----------------------------------------------------------------------------------------\n",t*defaultclock_dt);
			print_neurons(neurons, N_Group_T+N_S);
			SolveNeurons(neurons, N_Group_T, SpikeArray);	// maybe should bring the for inside out for(int i =0; i < N_T; i++) SolveNeuron(neurons[i],Spikearray[i]);
			printf("After SolveNeurons, t: %lf\n",t*defaultclock_dt);
			print_neurons(neurons, N_Group_T+N_S);
			/*printf("Loop %d\n",t);
			for(int i = 0; i < N_Group_T+N_S; i++){
				printf("%d\n",SpikeArray[i]);
			}
			*/
			//PoissonThreshold(input, N_S, N_Group_S, SpikeArray);
			if(t == 0 || t == 1) SpikeArray[25+N_Group_S] = 1;
			else SpikeArray[25+N_Group_S] = 0;
			
			if(t == 0) SpikeArray[46+N_Group_S] = 1;
			else SpikeArray[46+N_Group_S] = 0;

			//if(t*defaultclock_dt == 0.001) SpikeArray[25+N_Group_S] = 1;
			//else SpikeArray[25+N_Group_S] = 0;

			if(t == 2) SpikeArray[24+N_Group_S] = 1;
			else SpikeArray[24+N_Group_S] = 0;

			if(t == 2 || t == 13) SpikeArray[31+N_Group_S] = 1;
			else SpikeArray[31+N_Group_S] = 0;


			if(t == 5) SpikeArray[28+N_Group_S] = 1;
			else SpikeArray[28+N_Group_S] = 0;

			if(t == 6 || t == 7) SpikeArray[22+N_Group_S] = 1;
			else SpikeArray[22+N_Group_S] = 0;

			if(t == 8){
				SpikeArray[23+N_Group_S] = 1;
				SpikeArray[26+N_Group_S] = 1;
				SpikeArray[35+N_Group_S] = 1;
			} 
			else {
				SpikeArray[23+N_Group_S] = 0;
				SpikeArray[26+N_Group_S] = 0;
				SpikeArray[35+N_Group_S] = 0;
			}

			if(t == 9) SpikeArray[19+N_Group_S] = 1;
			else SpikeArray[19+N_Group_S] = 0;

			if(t == 11) SpikeArray[68+N_Group_S] = 1;
			else SpikeArray[68+N_Group_S] = 0;

			if(t == 12) SpikeArray[32+N_Group_S] = 1;
			else SpikeArray[32+N_Group_S] = 0;


			if(t == 14) SpikeArray[21+N_Group_S] = 1;
			else SpikeArray[21+N_Group_S] = 0;

			if(t == 14) SpikeArray[77+N_Group_S] = 1;
			else SpikeArray[77+N_Group_S] = 0;


			/*if(t*defaultclock_dt == 0.002) SpikeArray[6+N_Group_S] = 1;
			else SpikeArray[6+N_Group_S] = 0;

			if(t*defaultclock_dt == 0.002) SpikeArray[8+N_Group_S] = 1;
			else SpikeArray[8+N_Group_S] = 0;

			if(t*defaultclock_dt == 0.003) SpikeArray[1+N_Group_S] = 1;
			else SpikeArray[1+N_Group_S] = 0;*/

			printf("SpikeArray\n");
			for(int i = 0; i < N_Group_T+N_S; i++){
				printf("%d, ",SpikeArray[i]);
			}

			printf("\nSynapses//////////////////////////////////////////////////////\n");
			//print_synapses(syn,N_S+N_Group_S,N_Group_T+N_S);
			print_synapses(syn,N_S+N_Group_S,N_Group_S);
			UpdateSynapses_pre(syn, neurons, N_S+N_Group_S, N_Group_T+N_S, SpikeArray, t*defaultclock_dt);
			printf("\n\n\nSynapses after pre update, t: %lf\n\n\n",t*defaultclock_dt);
			
			//print_synapses(syn,N_S+N_Group_S,N_Group_T+N_S);
			print_synapses(syn,N_S+N_Group_S,N_Group_S);
			fflush(stdout);
			UpdateSynapses_post(syn, N_S+N_Group_S, N_Group_T, SpikeArray, t*defaultclock_dt);
			printf("\n\n\nSynapses after post update, t: %lf\n\n\n",t*defaultclock_dt);
			//print_synapses(syn,N_S+N_Group_S,N_Group_T+N_S);
			print_synapses(syn,N_S+N_Group_S,N_Group_S);
			print_neurons(neurons, N_Group_T+N_S);
		}
	}
	return 0;
}
