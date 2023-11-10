# M. Nicoletti et al. 2023, "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"

# This codes simulates RIM voltage and current clamp responses shown in Figure 7 panels A,D and G

# to simulate the knock-out set equal to zero the corresponding conductance


from neuron import h,gui
import numpy
from numpy import loadtxt
from matplotlib import pyplot

 
from VA5_simulation_iclamp import VA5_simulation_iclamp
from VA5_simulation_vclamp import VA5_simulation_vc

from g_to_Scm2 import gScm2

surf=389.3e-8 

#path='E:\RIM_MARCH_2021\VA5_MANUAL_MODEL_FINAL'

v=numpy.linspace(start=-60, stop=70, num=14)
ic=numpy.linspace(start=-0.01,stop=0.01,num=9)

names=numpy.array("slo2unc2,slo2egl19,slo2iso,EGL19,UNC2,cca1,irk,shk1mix,nca,leak,eleak,c2,cm".split(','))


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

ind=numpy.where(numpy.logical_and(best_time2[0]>=50, best_time2[0]<=60))
ind_max=numpy.amax(ind)
ind_min=numpy.amin(ind)
rp=numpy.mean(best_voltage[:,ind_min:ind_max],axis=1)
print('resting potential')
print(rp)



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
fig.legend(loc=5)
pyplot.title('IV-CURVES')
pyplot.show()













