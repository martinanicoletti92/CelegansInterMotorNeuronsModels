# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105


def VB6_simulation_iclamp(gVB6_scaled,s1,s2,ns):
    
 
    from neuron import h,gui
    import numpy
    import math


    surf=524.64e-8 # surface in cm^2 form neuromorpho VB6L
    vol=58.54e-12 # total volume
    L=math.sqrt(surf/math.pi)
    rsoma=L*1e4
    cm_uFcm2=1


    soma=h.Section(name="soma")
    soma.L=rsoma
    soma.diam=rsoma
    soma.Ra=100
    soma.cm=gVB6_scaled[15]
    h.psection(sec=soma)
    


    soma.insert('slo2egl19')
    soma.insert('slo1egl19')
    soma.insert('slo2unc2')
    soma.insert('slo1unc2')
    soma.insert('slo2iso')
    soma.insert('slo1iso')  
    soma.insert('egl19')
    soma.insert('unc2')
    soma.insert('cca1')    
    soma.insert('irk')
    soma.insert('shk1')
    soma.insert('leak')
    soma.insert('nca')
    soma.insert('cadiff')

    

       
    for seg in soma:
                

        
        seg.slo1egl19.gbar=gVB6_scaled[0]
        seg.slo2egl19.gbar = gVB6_scaled[1]
        seg.slo1unc2.gbar=gVB6_scaled[2]
        seg.slo2unc2.gbar = gVB6_scaled[3]        
        seg.slo1iso.gbar = gVB6_scaled[4]        
        seg.slo2iso.gbar=gVB6_scaled[5]               
        seg.egl19.gbar=gVB6_scaled[6] 
        seg.unc2.gbar=gVB6_scaled[7]
        seg.cca1.gbar=gVB6_scaled[8]
        seg.irk.gbar=gVB6_scaled[9]
        seg.shk1.gbar=gVB6_scaled[10]        
        seg.nca.gbar=gVB6_scaled[11]        
        seg.leak.gbar=gVB6_scaled[12]        
        seg.leak.e=gVB6_scaled[13]
        seg.slo2iso.c2=gVB6_scaled[14]

        seg.eca=60
        seg.ek=-80

    
    stim=h.IClamp(soma(0.5))
    dir(stim)
    
    stim.delay=5000
    stim.dur=1000
    
    v_vec = h.Vector()   
    t_vec = h.Vector() 
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
            

           
    v=[]
    v=numpy.array(list(ref_v))
    time1=numpy.array(ref_t)
    

    resc_ind=numpy.where(time1[1,:]>=4900)
    resc_min=numpy.amin(resc_ind)
    resc_max=numpy.amax(resc_ind)
    v_normalized=v[:,resc_min:resc_max]
    time=time1[:,resc_min:resc_max]-4900
    
  
  
    return v_normalized, time 
    
        


    
