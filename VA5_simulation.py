# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105


# VA5 voltage and current clamp responses shown in Figure 7 panels A,D and G

import os
from neuron import h,gui
import numpy
from numpy import loadtxt
from matplotlib import pyplot
from VA5_simulation_iclamp import VA5_simulation_iclamp
from VA5_simulation_vclamp import VA5_simulation_vc
from g_to_Scm2 import gScm2

surf=389.3e-8 

os.mkdir('VA5_SIMULATION')
path='VA5_SIMULATION'

v=numpy.linspace(start=-60, stop=70, num=14)
ic=numpy.linspace(start=-0.01,stop=0.01,num=9)


# CONDUCTANCES: slo2egl19,slo2iso,EGL19,irk,shk1,nca,leak,eleak,c2,cm
g0=[3,3,0.15,1,0.1,0.01,0.1,-70,1,1.5]

gstart=gScm2(g0,surf,6)

best_results=VA5_simulation_vc(gstart,-60,70,14)
best_current=numpy.array(list(best_results[0]))
best_iv=numpy.array(list(best_results[3]))
best_time=numpy.array(list(best_results[1]))
best_iv_peak=numpy.array(list(best_results[2]))
best_cai=numpy.array(list(best_results[4]))


best_cc=VA5_simulation_iclamp(gstart,-0.03,0.030,7)
best_voltage=best_cc[0]
best_time2=best_cc[1]



fig3=pyplot.figure(figsize=(8,4))
for i in range(0,14):
 curr_plot=pyplot.plot(best_time[i],best_current[i],color='red')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('I [pA]')
pyplot.title('Voltage clamp')
pyplot.show()


fig4=pyplot.figure(figsize=(8,4))
for i in range(0,7):
      volt_plot=pyplot.plot(best_time2[i],best_voltage[i],color='red')      
pyplot.xlabel('Time [ms]')
pyplot.ylabel('V [mV]')
pyplot.title('current_clamp')
pyplot.show()



fig=pyplot.figure(figsize=(8,4))
iv_plot=pyplot.plot(v,best_iv,color='red',marker='+',markersize=15)
pyplot.xlabel('mV')
pyplot.ylabel('pA')
pyplot.xlim(-130,100)
pyplot.title('IV-CURVES')
pyplot.show()













