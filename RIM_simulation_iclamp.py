# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105


def RIM_simulation_iclamp(gRIM_scaled,s1,s2,ns):
    
 
    from neuron import h,gui
    import numpy
    import math
    from operator import add

    
    
             
    surf=103.34e-8 # surface in cm^2 form neuromorpho RIML
    vol=10.94e-12 # total volume
    L=math.sqrt(surf/math.pi)
    rsoma=L*1e4
    cm_uFcm2=gRIM_scaled[8]


    soma=h.Section(name="soma")
    soma.L=rsoma
    soma.diam=rsoma
    soma.Ra=100
    soma.cm=cm_uFcm2
    h.psection(sec=soma)
    
    soma.insert('shl1')
    soma.insert('egl2')
    soma.insert('irk')
    soma.insert('cca1')
    soma.insert('unc2')
    soma.insert('egl19')
    soma.insert('leak')
    
        
    for seg in soma:
        
        seg.shl1.gbar=gRIM_scaled[0]
        seg.egl2.gbar=gRIM_scaled[1]
        seg.irk.gbar=gRIM_scaled[2]
        seg.cca1.gbar=gRIM_scaled[3]
        seg.unc2.gbar=gRIM_scaled[4]
        seg.egl19.gbar=gRIM_scaled[5]
        seg.leak.gbar=gRIM_scaled[6]
        
        seg.leak.e=gRIM_scaled[7]
           
        seg.eca=60
        seg.ek=-80
        
        
        
    stim=h.IClamp(soma(0.5))
    dir(stim)
    
    stim.delay=5000
    stim.amp=10
    stim.dur=5000
    
    v_vec = h.Vector()   
    t_vec = h.Vector()        # Time stamp vector
    v_vec.record(soma(0.5)._ref_v)
    t_vec.record(h._ref_t)

    simdur =14000

    ref_v=[]
    ref_t=[]

    num_step=11
    
    
    for i in numpy.linspace(start=s1, stop=s2, num=ns):
        
         stim.amp=i
         h.tstop=simdur
         h.dt=0.04
         h.finitialize(-60)
         h.run()
            
         ref_t_vec=numpy.zeros_like(t_vec)
         t_vec.to_python(ref_t_vec)
         ref_t.append(ref_t_vec)
            
        
         ref_v_vec=numpy.zeros_like(v_vec)
         v_vec.to_python(ref_v_vec)
         ref_v.append(ref_v_vec)
            
            # total current calculation
            
            
    v=[]
    v=numpy.array(list(ref_v))
    time1=numpy.array(ref_t)
    
    length=time1.shape
    
    # cut the initial transient 
    
    dd=numpy.amax(numpy.where(time1[1,:]<4000))
    
    time=time1[:,dd:length[1]]-4000

    volt=v[:, dd:length[1]]
    
        
    
    
    ## CALCULATION OF STEADY-STATE CURRENT-VOLATGE RELATION
    ind=numpy.where(numpy.logical_and(time[0]>=5990, time[0]<=6010))
    ind_max=numpy.amax(ind)
    ind_min=numpy.amin(ind)
    iv=numpy.mean(volt[:,ind_min:ind_max],axis=1)
	
	 # CALCULATION OF PEAK CURRENT-VOLTAGE RELATION (as in Ramot et al 2008)
    ind2=numpy.where(numpy.logical_and(time[0]>=1000, time[0]<=1100))
    ind2_max=numpy.amax(ind2)
    ind2_min=numpy.amin(ind2)
    iv_peak=numpy.amax(volt[:,ind2_min:ind2_max])
    iv_peak=[]
    

        
    for j in range(ns):
        if j<=3:
            peak=numpy.amin(volt[j,ind2_min:ind2_max])
        else:
            peak=numpy.amax(volt[j,ind2_min:ind2_max])
        iv_peak.append(peak)

    return volt, time, iv_peak, iv

    
    
        


