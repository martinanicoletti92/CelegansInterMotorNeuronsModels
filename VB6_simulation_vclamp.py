def VB6_simulation_vc(gVB6_scaled,vstart,vstop,ns):
    
# lista canali da cenGen tramite ricerca per geni 
# calcio e kcnl fisstao con i valori di t15.ode
    

    from neuron import h,gui
    import numpy
    import math

             
     
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
            
        
        
    stim=h.VClamp(soma(0.5))
    dir(stim)
    
    simdur = 1700
    stim.amp[0]=-60
    stim.amp[1]=-110
    stim.amp[2]=-60
    stim.dur[0]=1025
    stim.dur[1]=600
    stim.dur[2]=100
    
    ik_vec = h.Vector()   
    ica_vec=h.Vector()
    inca_vec=h.Vector()
    ileakVB6_vec=h.Vector()
    t_vec = h.Vector()  
    #curr_shl1_vec=h.Vector()

   
      # Time stamp vector
    ik_vec.record(soma(0.5)._ref_ik)
    ica_vec.record(soma(0.5)._ref_ica)
    inca_vec.record(soma(0.5)._ref_i_nca)
    ileakVB6_vec.record(soma(0.5)._ref_i_leak)
   # curr_shl1_vec.record(soma(0.5)._ref_curr_shl1)

    t_vec.record(h._ref_t)
    
    ref_ik=[]
    ref_ica=[]
    ref_t=[]
    ref_inca=[]
    ref_ileakVB6=[]
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
        ref_ileakVB6_vec=numpy.zeros_like(ileakVB6_vec)
        ileakVB6_vec.to_python(ref_ileakVB6_vec)
        ref_ileakVB6.append(ref_ileakVB6_vec)
    
    
    
     
    # total current calculation
    itot=[]
    #itot=map(sum, zip(ref_ik,ref_ica,ref_inca,ref_ileakVB6))   
    itot=map(sum, zip(ref_ik,ref_ica,ref_ileakVB6,ref_inca))        
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
