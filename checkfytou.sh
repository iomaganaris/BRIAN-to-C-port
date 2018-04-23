
make
echo "C running"
./simulation > neofytoutputC
echo "python running"
python prepostSTDP_savings.py > neofytoutputPython
cmp array_A.txt array_A_python.txt
cmp Neurons_I.txt Neuron_I_python.txt
cmp Spikes.txt spikespython.txt
