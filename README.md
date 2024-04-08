# CelegansInterMotorNeuronsModels

The folder contains the codes for reproducing the figures of the papers:
Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families
by M. Nicoletti et al. PloS ONE, 19(3): e0298105.
https://doi.org/10.1371/journal.pone.0298105

The repository contains python codes and NEURON .mod files for Hodgkin-Huxley models of
AVAL, AVAR, VA5, VB6, VD5, AIY, and RIM C.elegans neurons. 

The .mod files allow to model 22 ionic currents 
and the intracellular calcium. 

The pyhton files are organized for each neuron (AVAL, AVAR, AIY, RIM, VA5, VB6, VD5) as follows:
1. file for simulating voltage and current clamp responses entitled "NeuronName_simulation.py"
2.two files containing the functions for voltage and current clamp protocols entitled
"NeuronName_simulation_vclamp.py" and "NeuronName_simulation_iclamp.py", respectively

To reproduce the WT whole-cell behavior:
1. compile the NEURON .mod files in the same folder of the python codes
2. excute the file "NeuronName_simulation.py"

To simulate knock-out responses set to zero the conductance of the corresponding current. 
