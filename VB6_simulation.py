# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

# VB6 current and voltage clamp simulations shown in Figure 8 panels B,E, and H

from neuron import h,gui
import numpy
from numpy import loadtxt
from matplotlib import pyplot
import os
from VB6_simulation_iclamp import VB6_simulation_iclamp
from VB6_simulation_vclamp  import VB6_simulation_vc
from g_to_Scm2 import gScm2


os.mkdir('VB6_SIMULATION')
path='VB6_SIMULATION'

v=numpy.linspace(start=-60, stop=70, num=14)
ic=numpy.linspace(start=-0.01,stop=0.01,num=9)

surf=524.64e-8

# conductances: slo1EGL19,SLO2egl19,SLO2ISO,egl19,unc2,cca1,IRK,shk1mix,nca,leak,eleak,c2,cm
g0=[0.2,1.75,1,1,1,2,0.1,0.1,0,1,0.4,0.03,0.13,-52,1,1.5]

gstart=gScm2(g0,surf,12)

best_results=VB6_simulation_vc(gstart,-60,70,14)
best_current=numpy.array(list(best_results[0]))
best_iv=numpy.array(list(best_results[3]))
best_time=numpy.array(list(best_results[1]))
best_iv_peak=numpy.array(list(best_results[2]))


fname4="VB6_simulated_current_LEAK.txt"
fname5="VB6_simulated_timeVC_LEAK.txt"
fname6="VB6_simulated_IV_SS_LEAK.txt"

path4=os.path.join(path, fname4) 
path5=os.path.join(path, fname5) 
path6=os.path.join(path, fname6) 


numpy.savetxt(path4, best_current, delimiter="," , fmt="%s")
numpy.savetxt(path5, best_time, delimiter=", " , fmt="%s")
numpy.savetxt(path6, best_iv, delimiter="," , fmt="%s")



best_cc=VB6_simulation_iclamp(gstart,-0.03,0.030,7)
best_voltage=best_cc[0]
best_time2=best_cc[1]


ind=numpy.where(numpy.logical_and(best_time2[0]>=50, best_time2[0]<=60))
ind_max=numpy.amax(ind)
ind_min=numpy.amin(ind)
rp=numpy.mean(best_voltage[:,ind_min:ind_max],axis=1)
print('resting potential')
print(rp)


fname8='VB6_simulated_VOLTAGE_LEAK.txt'
fname9='VB6_simulated_time_CC_LEAK.txt'

path8=os.path.join(path, fname8) 
path9=os.path.join(path, fname9) 

numpy.savetxt(path8, best_voltage, delimiter="," , fmt="%s")
numpy.savetxt(path9, best_time2, delimiter=", " , fmt="%s")



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
pyplot.title('Current clamp')
pyplot.show()


fig=pyplot.figure(figsize=(8,4))
iv_plot=pyplot.plot(v,best_iv,color='red',marker='+',markersize=15)
pyplot.xlabel('V [mV]')
pyplot.ylabel('I [pA]')
pyplot.xlim(-130,80)
pyplot.title('IV STEADY STATE')
pyplot.show()


 












