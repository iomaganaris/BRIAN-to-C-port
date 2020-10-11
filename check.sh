cd build 
echo "C running"
./simulation > output
cd .. 
echo "python running"
python prepostSTDP_savings.py > output
cmp build/array_A.txt array_A_python.txt
cmp build/Neurons_I.txt Neuron_I_python.txt
cmp build/Spikes.txt spikespython.txt
