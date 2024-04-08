# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105


def VD5_simulation_vc(gVD5_scaled,vstart,vstop,ns):
    
# lista canali da cenGen tramite ricerca per geni 
# calcio e kcnl fisstao con i valori di t15.ode
    

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
        
        
    stim=h.VClamp(soma(0.5))
    dir(stim)
    
    simdur = 2400
    stim.amp[0]=-60
    stim.amp[1]=-110
    stim.amp[2]=-60
    stim.dur[0]=1025
    stim.dur[1]=1200
    stim.dur[2]=100
    
    ik_vec = h.Vector()   
    ica_vec=h.Vector()
    inca_vec=h.Vector()
    ileak_vec=h.Vector()
    t_vec = h.Vector()  
    #curr_shl1_vec=h.Vector()

   
      # Time stamp vector
    ik_vec.record(soma(0.5)._ref_ik)
    ica_vec.record(soma(0.5)._ref_ica)
    inca_vec.record(soma(0.5)._ref_i_nca)
    ileak_vec.record(soma(0.5)._ref_i_leak)
   # curr_shl1_vec.record(soma(0.5)._ref_curr_shl1)

    t_vec.record(h._ref_t)
    
    ref_ik=[]
    ref_ica=[]
    ref_t=[]
    ref_inca=[]
    ref_ileak=[]
    #ref_curr_shl1=[]

        
    for i in numpy.linspace(start=vstart, stop=vstop, num=ns):
        
        stim.amp[1]=i
        h.tstop=simdur
        h.dt=0.01
        h.finitialize(-60)
        h.run()
        
        #time
        ref_t_vec=numpy.zeros_like(t_vec)
        t_vec.to_python(ref_t_vec)
        ref_t.append(ref_t_vec)
        
        # potassium current
        ref_ik_vec=numpy.zeros_like(ik_vec)
        ik_vec.to_python(ref_ik_vec)
        ref_ik.append(ref_ik_vec)
        
              
        # shl1
        #ref_curr_shl1_vec=numpy.zeros_like(curr_shl1_vec)
        #curr_shl1_vec.to_python(ref_curr_shl1_vec)
        #ref_curr_shl1.append(ref_curr_shl1_vec)
        
               
        #calcium currents
        ref_ica_vec=numpy.zeros_like(ica_vec)
        ica_vec.to_python(ref_ica_vec)
        ref_ica.append(ref_ica_vec)
        
              
        # NCA currents
        ref_inca_vec=numpy.zeros_like(inca_vec)
        inca_vec.to_python(ref_inca_vec)
        ref_inca.append(ref_inca_vec)
    
        # LEAKAGE current
        ref_ileak_vec=numpy.zeros_like(ileak_vec)
        ileak_vec.to_python(ref_ileak_vec)
        ref_ileak.append(ref_ileak_vec)
    
    
    
     
    # total current calculation
    itot=[]
    #itot=map(sum, zip(ref_ik,ref_ica,ref_inca,ref_ileak))   
    itot=map(sum, zip(ref_ik,ref_ica,ref_ileak,ref_inca))        
    current=numpy.array(list(itot))
    inorm=current*1e9*surf #total current in pA
    #time array
    time1=numpy.array(ref_t)
    resc_ind=numpy.where(time1[1,:]>=1000)
    resc_min=numpy.amin(resc_ind)
    resc_max=numpy.amax(resc_ind)
    itot_normalized=inorm[:,resc_min:resc_max]
    time=time1[:,resc_min:resc_max]-1000
    
    
    ## CALCULATION OF STEADY-STATE CURRENT-VOLATGE RELATION
    ind=numpy.where(numpy.logical_and(time[0]>=600, time[0]<=625))
    ind_max=numpy.amax(ind)
    ind_min=numpy.amin(ind)
    iv=numpy.mean(itot_normalized[:,ind_min:ind_max],axis=1)
	
	 # CALCULATION OF PEAK CURRENT-VOLTAGE RELATION (as in Ramot et al 2008)
    ind2=numpy.where(numpy.logical_and(time[0]>=25, time[0]<=100))
    ind2_max=numpy.amax(ind2)
    ind2_min=numpy.amin(ind2)
    iv_peak=numpy.amax(itot_normalized[:,ind2_min:ind2_max])
    iv_peak=[]
    
    #for j in range(ns):
     #   peak=numpy.amax(itot_normalized[j,ind2_min:ind2_max])
     #   iv_peak.append(peak)
    #iv_peak=numpy.array(list(iv_peak))
    
    
    
    for j in range(ns):
        if j<=6:
            peak=numpy.amin(itot_normalized[j,ind2_min:ind2_max])
        else:
            peak=numpy.amax(itot_normalized[j,ind2_min:ind2_max])
        iv_peak.append(peak)

    return itot_normalized, time, iv_peak, iv
