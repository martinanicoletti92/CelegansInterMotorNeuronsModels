# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

def VD5_simulation_iclamp(gVD5_scaled,s1,s2,ns):
    
 
    from neuron import h,gui
    import numpy
    import math


    surf=351.53e-8 # surface in cm^2 form neuromorpho VD5L
    vol=51.3e-12 # total volume
    L=math.sqrt(surf/math.pi)
    rsoma=L*1e4
    cm_uFcm2=1


    soma=h.Section(name="soma")
    soma.L=rsoma
    soma.diam=rsoma
    soma.Ra=100
    soma.cm=gVD5_scaled[10]
    h.psection(sec=soma)
    



    soma.insert('slo2egl19')
    soma.insert('slo2iso')
    soma.insert('egl19')
    soma.insert('cca1')
    soma.insert('irk')
    soma.insert('shk1')
    soma.insert('leak')
    soma.insert('nca')
    soma.insert('cadiff')

    

       
    for seg in soma:
                

        
        seg.slo2egl19.gbar = gVD5_scaled[0]        
        seg.slo2iso.gbar=gVD5_scaled[1] 
        seg.egl19.gbar=gVD5_scaled[2] 
        seg.cca1.gbar=gVD5_scaled[3]
        seg.irk.gbar=gVD5_scaled[4]
        seg.shk1.gbar=gVD5_scaled[5]        
        seg.nca.gbar=gVD5_scaled[6]        
        seg.leak.gbar=gVD5_scaled[7]        
        seg.leak.e=gVD5_scaled[8]
        seg.slo2iso.c2=gVD5_scaled[9]


        
        seg.eca=60
        seg.ek=-80

    
    stim=h.IClamp(soma(0.5))
    dir(stim)
    
    stim.delay=5000
    stim.amp=10
    stim.dur=1000
    
    v_vec = h.Vector()   
    t_vec = h.Vector()        # Time stamp vector
    v_vec.record(soma(0.5)._ref_v)
    t_vec.record(h._ref_t)

    simdur =7000

    ref_v=[]
    ref_t=[]

    
    
    for i in numpy.linspace(start=s1, stop=s2, num=ns):
        
         stim.amp=i
         h.tstop=simdur
         h.dt=0.4
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
    

    resc_ind=numpy.where(time1[1,:]>=4900)
    resc_min=numpy.amin(resc_ind)
    resc_max=numpy.amax(resc_ind)
    v_normalized=v[:,resc_min:resc_max]
    time=time1[:,resc_min:resc_max]-4900
    
    
    ## CALCULATION OF STEADY-STATE CURRENT-VOLATGE RELATION
#     ind=numpy.where(numpy.logical_and(time[0]>=5060, time[0]<=5100))
#     ind_max=numpy.amax(ind)
#     ind_min=numpy.amin(ind)
#     iv=numpy.mean(v_normalized[:,ind_min:ind_max],axis=1)
# 	
# 	 # CALCULATION OF PEAK CURRENT-VOLTAGE RELATION (as in Ramot et al 2008)
#     ind2=numpy.where(numpy.logical_and(time[0]>=100, time[0]<=300))
#     ind2_max=numpy.amax(ind2)
#     ind2_min=numpy.amin(ind2)
#     iv_peak=numpy.amax(v_normalized[:,ind2_min:ind2_max])
#     iv_peak=[]
    

    
    
    # for j in range(ns):
    #     if j<=2:
    #         peak=numpy.amin(v_normalized[j,ind2_min:ind2_max])
    #     else:
    #         peak=numpy.amax(v_normalized[j,ind2_min:ind2_max])
    #     iv_peak.append(peak)

    return v_normalized, time 
    
        


    
