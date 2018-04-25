cd /home/ioannis/Desktop/Porting
make
echo "C running"
./simulation > output
cd /home/ioannis/Desktop/Diploma_Thesis/Tests/Debbuging_code
echo "python running"
python prepostSTDP_savings.py > output
cmp /home/ioannis/Desktop/Porting/array_A.txt /home/ioannis/Desktop/Diploma_Thesis/Tests/Debbuging_code/array_A_python.txt
cmp /home/ioannis/Desktop/Porting/Neurons_I.txt /home/ioannis/Desktop/Diploma_Thesis/Tests/Debbuging_code/Neuron_I_python.txt
cmp /home/ioannis/Desktop/Porting/Spikes.txt /home/ioannis/Desktop/Diploma_Thesis/Tests/Debbuging_code/spikespython.txt
