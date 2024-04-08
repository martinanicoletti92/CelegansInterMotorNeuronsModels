# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105


# AVAL neuron H-H MODEL 
# current and voltage clamp simulations shown in Figure 1 panels A, C, D, and F

import numpy
from numpy import loadtxt
from matplotlib import pyplot
import os
from AVAL_simulation_iclamp import AVA_simulation_iclamp
from AVAL_simulation_vclamp import AVA_simulation_vc
from g_to_Scm2 import gScm2

os.mkdir('AVAL_SIMULATION')
path='AVAL_SIMULATION'


v=numpy.linspace(start=-110, stop=50, num=17)
surf=1123.84e-8# 


# coductances: egl19, leak, irk, nca, eleak, cm
g0=[0.104385,0.150164,0.1,0,-39,0.859551]

gbest=gScm2(g0,surf,3)


best_cc=AVA_simulation_iclamp(gbest,-0.03,0.03,7)
best_voltage=best_cc[0]
best_time2=best_cc[1]
best_vipeaks=best_cc[2]
best_viss=best_cc[3]


fname8="AVAL_simulated_VOLTAGE_WT.txt"
fname9="AVAL_simulated_timeCC_WT.txt"
fname10="AVAL_simulated_VI_SS_WT.txt"
fname11="AVAL_IV_simulated_VI_PEAKS_WT.txt"


path8=os.path.join(path, fname8) 
path9=os.path.join(path, fname9) 
path10=os.path.join(path, fname10) 
path11=os.path.join(path, fname11) 


numpy.savetxt(path8, best_voltage, delimiter="," , fmt="%s")
numpy.savetxt(path9, best_time2, delimiter=", " , fmt="%s")
numpy.savetxt(path10, best_viss, delimiter="," , fmt="%s")
numpy.savetxt(path11, best_vipeaks, delimiter=", " , fmt="%s")


best_results=AVA_simulation_vc(gbest,-110,50,17)
best_current=numpy.array(list(best_results[0]))
best_iv=numpy.array(list(best_results[3]))
best_time=numpy.array(list(best_results[1]))
best_iv_peak=numpy.array(list(best_results[2]))



fname4="AVAL_simulated_current_WT.txt"
fname5="AVAL_simulated_timeVC_WT.txt"
fname6="AVAL_simulated_IV_SS_WT.txt"
fname7="AVAL_IV_simulated_PEAKS_WT.txt"



path4=os.path.join(path, fname4) 
path5=os.path.join(path, fname5) 
path6=os.path.join(path, fname6) 
path7=os.path.join(path, fname7) 

numpy.savetxt(path4, best_current, delimiter="," , fmt="%s")
numpy.savetxt(path5, best_time, delimiter="," , fmt="%s")
numpy.savetxt(path6, best_iv, delimiter=", " , fmt="%s")
numpy.savetxt(path7, best_iv_peak, delimiter=", " , fmt="%s")


fig4=pyplot.figure(figsize=(8,4))
for i in range(0,7):
 volt_plot=pyplot.plot(best_time2[i],best_voltage[i],color='red',linestyle='solid')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('V [mV]')
pyplot.title('Current_Clamp')
pyplot.show()




fig3=pyplot.figure(figsize=(8,4))
for i in range(0,17):
 curr_plot=pyplot.plot(best_time[i],best_current[i],color='red',linestyle='solid')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('I [pA]')
pyplot.title('Voltage clamp')
pyplot.show()




# plot

fig=pyplot.figure(figsize=(8,4))
iv_plot=pyplot.plot(v,best_iv,color='red',marker='+',markersize=15)
pyplot.xlabel('V [mV]')
pyplot.ylabel('I [pA]')
pyplot.xlim(-130,60)
pyplot.title('IV steady-state')
pyplot.show()

  


fig2=pyplot.figure(figsize=(8,4))
iv2_plot=pyplot.plot(v,best_iv_peak,color='red',marker='+',markersize=15)
pyplot.xlabel('V [mV]')
pyplot.ylabel('I [pA]')
pyplot.xlim(-130,60)
pyplot.title('IV PEAKS')
pyplot.show()
