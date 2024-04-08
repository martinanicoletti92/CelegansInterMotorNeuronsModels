# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

# AIY neuron H-H MODEL
# current and voltage clamp simulations shown in Figure 3

import os
from neuron import h,gui
import numpy
from matplotlib import pyplot
from AIY_simulation_iclamp import AIY_simulation_iclamp
from AIY_simulation_vclamp import AIY_simulation_vc
from g_to_Scm2 import gScm2


os.mkdir('AIY_SIMULATION')
path='AIY_SIMULATION'

v=numpy.linspace(start=-120, stop=50, num=18)
ic=numpy.linspace(start=-15, stop=35, num=11)
surf=65.89e-8# surface in cm^2 form neuromorpho AIYL

#conductances: leak, slo1iso,kqt1,egl19,slo1egl19,nca,irk,eleak,cm
g0=[0.14,1,0.2,0.1,0.92,0.06,0.5,-89.57,1.6]


gbest=gScm2(g0,surf,6)


best_results=AIY_simulation_vc(gbest,-120,50,18)
best_current=numpy.array(list(best_results[0]))
best_iv=numpy.array(list(best_results[3]))
best_time=numpy.array(list(best_results[1]))
best_iv_WT=numpy.array(list(best_results[2]))

fname4="AIY_simulated_current_WT.txt"
fname5="AIY_simulated_time_WT.txt"
fname6="AIY_simulated_IV_SS_WT.txt"
fname7="AIY_IV_simulated_PEAK_WT.txt"


path4=os.path.join(path, fname4) 
path5=os.path.join(path, fname5) 
path6=os.path.join(path, fname6) 
path7=os.path.join(path, fname7) 

numpy.savetxt(path4, best_current, delimiter="," , fmt="%s")
numpy.savetxt(path5, best_time, delimiter="," , fmt="%s")
numpy.savetxt(path6, best_iv, delimiter=", " , fmt="%s")
numpy.savetxt(path7, best_iv_WT, delimiter=", " , fmt="%s")


best_cc=AIY_simulation_iclamp(gbest,-0.015,0.035,11)
best_voltage=best_cc[0]
best_time2=best_cc[1]
best_VIss=best_cc[3]
best_VIpeaks=best_cc[2]

fname8='AIY_CC_simulated_voltage_WT.txt'
fname9='AIY_CC_simulated_time_WT.txt'
path8=os.path.join(path, fname8) 
path9=os.path.join(path, fname9) 

numpy.savetxt(path8, best_voltage, delimiter="," , fmt="%s")
numpy.savetxt(path9, best_time2, delimiter=", " , fmt="%s")

# plot

fig=pyplot.figure(figsize=(8,4))
iv_plot=pyplot.plot(v,best_iv,color='red',marker='+',markersize=15,label='optimized-ss')
iv_plot=pyplot.plot(v,best_iv_WT,color='red',marker='o',markersize=15,label='optimized-peaks')
pyplot.xlabel('V [mV]')
pyplot.ylabel('I [pA]')
pyplot.xlim(-130,60)
fig.legend(loc=5)
pyplot.title('IV-CURVES')
pyplot.show()



fig3=pyplot.figure(figsize=(8,4))
for i in range(0,18):
 curr_plot=pyplot.plot(best_time[i],best_current[i],color='red',linestyle='solid')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('I [pA]')
pyplot.title('Voltage clamp')
pyplot.show()




fig4=pyplot.figure(figsize=(8,4))
for i in range(0,10):
 volt_plot=pyplot.plot(best_time2[i],best_voltage[i],color='red',linestyle='solid')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('V [mV]')
#pyplot.xlim(0,0.7)
pyplot.title('Current_Clamp')
pyplot.show()
             

