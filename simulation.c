/**
 * @file simulation.c
 * @author Ioannis Magkanaris
 * @author Alexandros Neofytou
 * @date 23 April 2018
 * @brief Simulation main file.
 *
 * Contains the main function of the program, with all the initializations
 * of the variables needed for a specific simulation.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//#include "neuron.h"
#include "synapses.h"


#define NxM
//#define MxM
//#define NxMxM	//not working probably
//#define allconnected

double defaultclock_dt = 1*1e-3;	//ms

//double taum_S = 10 * 1e-3; //ms
double Ee = 0 * 1e-3; //mV
double tuae = 2 * 1e-3; //ms
double Fon = 50; //Hz
double Foff = 3; //Hz

#ifdef NxM 
	double s = 100 * 1e-10;//100*1e-10; 
#endif
#ifdef MxM 
	double s = 500000;	//for testing
#endif
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

int main(void){

	FILE *f = fopen("Spikes.txt", "w");
	FILE *g = fopen("Neurons_I.txt", "w");
	FILE *h = fopen("array_A.txt", "w");
//neofytou connectivity testing
	FILE *in = fopen("connections.txt", "r");
	int nruns = 1;

	for(int nrun = 0; nrun < nruns; nrun++){
		double realtime = 0;
		double stime = 0.1; //second
		double stime2 = 50; //second

		double resolution_export = 10 * 1e-3; //every x ms

		int N = 100;
		#ifdef MxM
			N = 0;
		#endif
		int M = 2;

		int N_S;//100;
		int N_Group_S;
		int N_Group_T;

		#ifdef NxM
			N_S = N;//100;
			N_Group_S = 0;
			N_Group_T = M;	//for the simulation we have, normally is N
	    	#endif

		#ifdef MxM
			N_S = 0;//100;
			N_Group_S = M;
			N_Group_T = M;	//for the simulation we have, normally is N
		#endif
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
	    //printf("F_input1\n");
	    for(int i = 0; i < N_S; i++){
	    	F_input1[i] = Foff;			//maybe is not needed
	    	F_input1[i] = exp(-(pow((i+1)-input1_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	    	//printf("%lf ",F_input1[i]);
	    }
	    //printf("\n");

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
	    neurons = (Neuron*)malloc(sizeof(Neuron)*(N_Group_T));
	    for(int i = 0; i<N_Group_T; i++){				// Initiliazation of Neurons
	    	neurons[i].vt = vtrest;
	    	#ifdef NxM
	    		neurons[i].vm = EL;
	    	#endif
	    	#ifdef MxM
				neurons[i].vm = vtrest + 0.005;//EL;
			#endif
	    	//if(i%2==1) neurons[i].vm = vtrest + 0.005;
	    	neurons[i].I = 0;
	    	neurons[i].x = 0;
	    	neurons[i].Spike = 0;
	    }

	    int *SpikeArray;
	    SpikeArray = (int*)malloc(sizeof(int)*(N_Group_T+N_S));
	    
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
	    for(int i = 0; i < N_S+N_Group_S; i++) syn[i] = (Synapse*)malloc(sizeof(Synapse) * (N_Group_T));

//Connectivity
	    #ifndef allconnected
	    int con = 0;
	       
	    for(int i = 0; i < N_Group_S; i++){
	    	for(int j = 0; j < N_Group_T; j++){ 
				fscanf (in, "%d", &con);
	      	   	if (con == 1) syn[i][j].conn = 1;
				else syn[i][j].conn = 0; 
			}    
    	}

	    for(int i = N_Group_S; i < N_Group_S+N_S; i++){
	    	for(int j = 0; j < N_Group_T; j++){ 
				fscanf (in, "%d", &con);
	      	   	if (con == 1) syn[i][j].conn = 1;
				else syn[i][j].conn = 0; 
			}    
    	}
    	#else
    	for(int i = 0; i < N_Group_S+N_S; i++)
	    	for(int j = 0; j < N_Group_T; j++)
	    		syn[i][j].conn = 1;
		#endif

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
    	int init_const = 0;
	    for(int i = 0; i < N_Group_S; i++){
	    	for(int j = 0; j < N_Group_T; j++){
			    if (syn[i][j].conn) {
		    		//syn[i][j].conn = 1;	// all connected
	//Connectivity, initialization now happens only in connected synapses.
		    		syn[i][j].FBp = 0;
		    		syn[i][j].FBn = 0;
		    		syn[i][j].R = 1;
		    		#ifdef MxM
		    			syn[i][j].u = 1;	// for testing
		    		#endif
		    		syn[i][j].U = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
		    		syn[i][j].A = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
		    		init_const++;
		    		//syn[i][j].U = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
		    		//syn[i][j].A = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
			    }
	    	}
	    }
	    // Initialization of Synapses for external input (bottom rows)
	    for(int i = N_Group_S; i < N_Group_S+N_S; i++){
	        for(int j = 0; j < N_Group_T; j++){
			    if (syn[i][j].conn) {
	//neofytou connectivity
			    	//syn[i][j].conn = 1;	// all connected
		    		syn[i][j].FBp = 0;
		    		syn[i][j].FBn = 0;
		    		syn[i][j].R = 1;
		    		//syn[i][j].U = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
		    		//syn[i][j].A = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
		    		syn[i][j].U = exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
		    		syn[i][j].A = exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
			     	init_const++;
			    }
	        }
	     }
	    // Initialization of Synapses for external input ( right columns)
	    /*for(int i = 0; i < N_S + N_Group_S; i++){
		    for(int j = N_Group_T; j < N_Group_T+N_S; j++){
		    	syn[i][j].conn = 1;	// all connected
	    		syn[i][j].FBp = 0;
	    		syn[i][j].FBn = 0;
	    		syn[i][j].R = 1;
	    		syn[i][j].U = exp(-(((pow(((i*N_Group_T+j-N_Group_T)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
	    		syn[i][j].A = exp(-(((pow(((i*N_Group_T+j-N_Group_T)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
		    }
		}*/

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
			//printf("t: %.20lf----------------------------------------------------------------------------------------\n",t*defaultclock_dt);
			//print_neurons(neurons, N_Group_T);
			fprintf(g, "t: %.3lf\n", t*defaultclock_dt);
			for(int i = 0; i < N_Group_T; i++){
				fprintf(g,"%.8e ",neurons[i].I);
			}
			fprintf(g, "\n");
			SolveNeurons(neurons, N_Group_T, SpikeArray);	// maybe should bring the for inside out for(int i =0; i < N_T; i++) SolveNeuron(neurons[i],Spikearray[i]);
			//printf("After SolveNeurons, t: %lf\n",t*defaultclock_dt);
			//print_neurons(neurons, N_Group_T);
			/*printf("Loop %d\n",t);
			for(int i = 0; i < N_Group_T+N_S; i++){
				printf("%d\n",SpikeArray[i]);
			}
			*/
			//PoissonThreshold(input, N_S, N_Group_S, SpikeArray);
			#ifdef NxM
				//if(t == 0 ) SpikeArray[89+N_Group_T] = 1;
				//else SpikeArray[89+N_Group_T] = 0;
				
				//if(t == 2 || t == 25) SpikeArray[73+N_Group_T] = 1;
				//else SpikeArray[73+N_Group_T] = 0;

				//if(t*defaultclock_dt == 0.001) SpikeArray[25+N_Group_S] = 1;
				//else SpikeArray[25+N_Group_S] = 0;

				if(t == 4) SpikeArray[23+N_Group_T] = 1;
				else SpikeArray[23+N_Group_T] = 0;
				
				if(t == 5){
					SpikeArray[20+N_Group_T] = 1;
					SpikeArray[45+N_Group_T] = 1;
					
				} 
				else {
					SpikeArray[20+N_Group_T] = 0;
					SpikeArray[45+N_Group_T] = 0;
					
				}

				//if(t == 6 ) SpikeArray[81+N_Group_T] = 1;
				//else SpikeArray[81+N_Group_T] = 0;
				
				if(t == 7 ) SpikeArray[7+N_Group_T] = 1;
				else SpikeArray[7+N_Group_T] = 0;
				
				if(t == 7 || t == 11 || t == 13 || t == 20 || t == 26 || t == 28 || t == 31) SpikeArray[27+N_Group_T] = 1;
				else SpikeArray[27+N_Group_T] = 0;

				if(t == 9 ) SpikeArray[25+N_Group_T] = 1;
				else SpikeArray[25+N_Group_T] = 0;
				
				if(t == 11 ) SpikeArray[11+N_Group_T] = 1;
				else SpikeArray[11+N_Group_T] = 0;
				
				if(t == 15 ) SpikeArray[29+N_Group_T] = 1;
				else SpikeArray[29+N_Group_T] = 0;
				
				if(t == 19 ) SpikeArray[32+N_Group_T] = 1;
				else SpikeArray[32+N_Group_T] = 0;
				
				if(t == 25){
					SpikeArray[41+N_Group_T] = 1;
					//SpikeArray[73+N_Group_T] = 1;
					
				} 
				else {
					SpikeArray[41+N_Group_T] = 0;
					//SpikeArray[73+N_Group_T] = 0;
					
				}
				
				if(t == 26 ) SpikeArray[14+N_Group_T] = 1;
				else SpikeArray[14+N_Group_T] = 0;
				
				if(t == 28 || t == 30 ) SpikeArray[22+N_Group_T] = 1;
				else SpikeArray[22+N_Group_T] = 0;
				
				if(t == 28 ) SpikeArray[31+N_Group_T] = 1;
				else SpikeArray[31+N_Group_T] = 0;
				
				if(t == 29 ) SpikeArray[24+N_Group_T] = 1;
				else SpikeArray[24+N_Group_T] = 0;
				
			#endif


			/*if(t*defaultclock_dt == 0.002) SpikeArray[6+N_Group_S] = 1;
			else SpikeArray[6+N_Group_S] = 0;
			if(t*defaultclock_dt == 0.002) SpikeArray[8+N_Group_S] = 1;
			else SpikeArray[8+N_Group_S] = 0;
			if(t*defaultclock_dt == 0.003) SpikeArray[1+N_Group_S] = 1;
			else SpikeArray[1+N_Group_S] = 0;*/
			int flag = 0;
			//printf("SpikeArray\n");
			for(int i = 0; i < N_Group_T+N_S; i++){
				if (SpikeArray[i]==1) flag=1;
				//printf("%d, ",SpikeArray[i]);
			}
			fprintf(f, "t: %.3lf \n", t*defaultclock_dt);
			for(int i = 0; i < N_Group_T; i++){
				if (SpikeArray[i]!=0){
					fprintf(f,"%d ",i);
				}
			}
			fprintf(f, "\n");
			if(N_S>0)fprintf(f, "t: %.3lf \n", t*defaultclock_dt);
			for(int i = N_Group_T; i < N_S; i++){
				if(SpikeArray[i]!=0){
					fprintf(f,"%d ",i-N_Group_T);
				}
			}
			if(N_S>0) fprintf(f, "\n");
			
			//if(N_S>0){
			
			//}
			printf("\nSynapses//////////////////////////////////////////////////////\n");
			//print_synapses(syn,N_S+N_Group_S,N_Group_T+N_S);
			//print_synapses(syn,N_S+N_Group_S,N_Group_T);
			flag = 0;
			for(int i = N_Group_T; i < N_Group_T+N_S; i++){
				if(SpikeArray[i]==1)flag=1;
			}
			if(flag==1){
				fprintf(h, "t: %.3lf \n", t*defaultclock_dt);
				for(int i =0; i < N_S+N_Group_S; i++){
					for(int j = 0; j < N_Group_T; j++){
//neofytou connectivity
						if (syn[i][j].conn) fprintf(h,"%.8e ", syn[i][j].A);	
						//if((i*N_T+j+1)%4 == 0) printf("\n");
					}
					//printf("\n");
				}
				fprintf(h, "\n");
			}
			UpdateSynapses_pre(syn, neurons, N_S, N_Group_S, N_Group_T, SpikeArray, t*defaultclock_dt);
			//printf("\n\n\nSynapses after pre update, t: %lf\n\n\n",t*defaultclock_dt);
			flag = 0;
			for(int i = 0; i < N_Group_T; i++){
				if(SpikeArray[i]==1)flag=1;
			}
			//print_synapses(syn,N_S+N_Group_S,N_Group_T+N_S);
			//print_synapses(syn,N_S+N_Group_S,N_Group_T);
			if(flag==1){
				fprintf(h, "t: %.3lf \n", t*defaultclock_dt);
				for(int i =0; i < N_S+N_Group_S; i++){
					for(int j = 0; j < N_Group_T; j++){
//neofytou connectivity
						if (syn[i][j].conn) fprintf(h,"%.8e ", syn[i][j].A);	
						//if((i*N_T+j+1)%4 == 0) printf("\n");
					}
					//printf("\n");
				}
				fprintf(h, "\n");
			}
			fflush(stdout);
			UpdateSynapses_post(syn, N_S, N_Group_S, N_Group_T, SpikeArray, t*defaultclock_dt);
			//printf("\n\n\nSynapses after post update, t: %lf\n\n\n",t*defaultclock_dt);
			//print_synapses(syn,N_S+N_Group_S,N_Group_T+N_S);
			//print_synapses(syn,N_S+N_Group_S,N_Group_T);
			//print_neurons(neurons, N_Group_T);
		}
	}
	fclose(f);
	fclose(g);
	fclose(h);
//neofytou connectivity 
	fclose(in);
	return 0;
}
