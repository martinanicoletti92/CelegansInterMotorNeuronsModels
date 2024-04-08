# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

# This codes simulates RIM voltage and current clamp responses shown in Figure 5

import os
import numpy
from neuron import h,gui
from numpy import loadtxt
from matplotlib import pyplot
from RIM_simulation_iclamp import RIM_simulation_iclamp
from RIM_simulation_vclamp import RIM_simulation_vc


os.mkdir('RIM_SIMULATION')
path='RIM_SIMULATION'

v=numpy.linspace(start=-100, stop=50, num=16)
ic=numpy.linspace(start=-15,stop=35,num=11)

surf=103.34e-8 # surface in cm^2 form neuromorpho RIML


# conductances in S/cm^2: SHL1, EGL2, IRK, CCA1, unc2, egl19, LEAK,eleak, cm

g=[0.0009048750067326097,0.0001411644285181245,0.0003272854640954744,0.0008451919806776876,9.676795045480941e-05,0.00032005818627638106,9.676795045480941e-05,-50,1.5]


# ============= SIMULATIONS ======================================

best_cc=RIM_simulation_iclamp(g,-0.015,0.035,11)
best_voltage=best_cc[0]
best_time2=best_cc[1]
best_vipeaks=best_cc[2]
best_viss=best_cc[3]


fname8='RIM_simulated_VOLTAGE.txt'
fname9='RIM_simulated_time_CC.txt'

path8=os.path.join(path, fname8) 
path9=os.path.join(path, fname9) 

numpy.savetxt(path8, best_voltage, delimiter="," , fmt="%s")
numpy.savetxt(path9, best_time2, delimiter=", " , fmt="%s")


best_results=RIM_simulation_vc(g,-100,50,16)
best_current=numpy.array(list(best_results[0]))
best_iv=numpy.array(list(best_results[3]))
best_time=numpy.array(list(best_results[1]))
best_iv_peak=numpy.array(list(best_results[2]))


fname4="RIM_simulated_current.txt"
fname5="RIM_simulated_time.txt"
fname6="RIM_simulated_IV_SS.txt"
fname7="RIM_IV_simulated_PEAKS.txt"

path4=os.path.join(path, fname4) 
path5=os.path.join(path, fname5) 
path6=os.path.join(path, fname6) 
path7=os.path.join(path, fname7) 

numpy.savetxt(path4, best_current, delimiter="," , fmt="%s")
numpy.savetxt(path5, best_time, delimiter=", " , fmt="%s")


# ================== FIGURES ===================================


fig=pyplot.figure(figsize=(8,4))
iv_plot=pyplot.plot(v,best_iv,color='red',marker='+',markersize=15,label='Model-SS')
iv3_plot=pyplot.plot(v,best_iv_peak,color='red',marker='+',markersize=15,label='Model-Peaks')
pyplot.xlabel('V [mV]')
pyplot.ylabel('I [pA]')
pyplot.xlim(-130,60)
fig.legend(loc=5)
pyplot.title('RIM IV-CURVES')
pyplot.show()



fig2=pyplot.figure(figsize=(8,4))
iv3_plot=pyplot.plot(ic,best_vipeaks,color='red',marker='+',markersize=15,label='Model-Peaks')
iv3_plot=pyplot.plot(ic,best_viss,color='red',marker='+',markersize=15,label='Model-SS')
pyplot.ylabel('V [mv]')
pyplot.xlabel('I [pA]')
pyplot.xlim(-130,60)
fig2.legend(loc=5)
pyplot.title('RIM VI-CURVES')



fig4=pyplot.figure(figsize=(8,4))
for i in range(0,10):
 volt_plot=pyplot.plot(best_time2[i],best_voltage[i],color='red',label='Model')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('V [mV]')
pyplot.title('RIM current-clamp')
pyplot.show()



fig3=pyplot.figure(figsize=(8,4))
for i in range(0,15):
 curr_plot=pyplot.plot(best_time[i],best_current[i],color='red',label='Model')
pyplot.xlabel('Time [ms]')
pyplot.ylabel('I [pA]')
pyplot.title('RIM voltage-clamp')
pyplot.show()





