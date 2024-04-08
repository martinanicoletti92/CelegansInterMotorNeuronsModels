# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

def AVA_simulation_iclamp(gAVA_scaled,s1,s2,ns):
    
    from neuron import h,gui
    import numpy
    import math
 
   
    
    surf=1123.84e-8 # surface in cm^2 form neuromorpho AVAL
    vol=129.6e-12 # total volume
    L=math.sqrt(surf/math.pi)
    rsoma=L*1e4

    cm_uFcm2=gAVA_scaled[5]
    cm_pF=cm_uFcm2*1e6/surf

   
    soma=h.Section(name="soma")
    soma.L=rsoma
    soma.diam=rsoma
    soma.Ra=100
    soma.cm=cm_uFcm2
    h.psection(sec=soma)

    soma.insert('irk')
    soma.insert('leak')
    soma.insert('egl19')
    soma.insert('nca')
   
    
    
    for seg in soma:
        

        seg.egl19.gbar=gAVA_scaled[0]
        seg.leak.gbar=gAVA_scaled[1]
        seg.irk.gbar=gAVA_scaled[2]
        seg.nca.gbar=gAVA_scaled[3]
        seg.leak.e=gAVA_scaled[4]
        
       
        seg.eca=60
        seg.ek=-80
        
        
    stim=h.IClamp(soma(0.5))
    dir(stim)
    
    stim.delay=1023
    stim.amp=10
    stim.dur=1000
    
    v_vec = h.Vector()   
    t_vec = h.Vector()       
    v_vec.record(soma(0.5)._ref_v)
    t_vec.record(h._ref_t)

    simdur =2500

    ref_v=[]
    ref_t=[]

    
    
    for i in numpy.linspace(start=s1, stop=s2, num=ns):
        
         stim.amp=i
         h.tstop=simdur
         h.dt=0.025
         h.finitialize(-60)
         h.run()
            
         ref_t_vec=numpy.zeros_like(t_vec)
         t_vec.to_python(ref_t_vec)
         ref_t.append(ref_t_vec)
            
        
         ref_v_vec=numpy.zeros_like(v_vec)
         v_vec.to_python(ref_v_vec)
         ref_v.append(ref_v_vec)
            

            
            
    v=[]
    v=numpy.array(list(ref_v))
    time1=numpy.array(ref_t)
    

    resc_ind=numpy.where(time1[1,:]>=1000)
    resc_min=numpy.amin(resc_ind)
    resc_max=numpy.amax(resc_ind)
    v_normalized=v[:,resc_min:resc_max]
    time=time1[:,resc_min:resc_max]-1000
    
    
    ## SS V-I curve
    ind=numpy.where(numpy.logical_and(time[0]>=23, time[0]<=63))
    ind_max=numpy.amax(ind)
    ind_min=numpy.amin(ind)
    iv=numpy.mean(v_normalized[:,ind_min:ind_max],axis=1)
	
	# PEAKS V-I curve
    ind2=numpy.where(numpy.logical_and(time[0]>=953, time[0]<=1023))
    ind2_max=numpy.amax(ind2)
    ind2_min=numpy.amin(ind2)
    iv_peak=numpy.amax(v_normalized[:,ind2_min:ind2_max])
    iv_peak=[]
    
  
    
    for j in range(ns):
        if j<=3:
            peak=numpy.amin(v_normalized[j,ind2_min:ind2_max])
        else:
            peak=numpy.amax(v_normalized[j,ind2_min:ind2_max])
        iv_peak.append(peak)

    return v_normalized, time, iv_peak, iv    
    
    


    
    
        


