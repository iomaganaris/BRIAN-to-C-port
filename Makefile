simulation_make: simulation.c neuron.c neuron.h synapses.c synapses.h
	gcc simulation.c neuron.c synapses.c -lm -o simulation
