
# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

def AIY_simulation_vc(gAIY_scaled,vstart,vstop,ns):
    
    # voltage clamp per Minimizzazioni con canali calcio fissati e ottimizzazione sul potassio
    #i valori delle conduttanze del calcio vengono dal set xpp DEFINITIVO.ode
    #fissati anche i valori di irk, kcnl per fittare le code e di slo1 slo2
    #libere anche le conduttanze di perdita
    
   
    from neuron import h,gui
    import numpy
    import math
    from matplotlib import pyplot
    from operator import add
    
         
    from neuron import h,gui
    import numpy
    import math
    from matplotlib import pyplot
    from operator import add

    surf=65.89e-8 # surface in cm^2 form neuromorpho AIYL
    vol=7.42e-12 # total volume
    L=math.sqrt(surf/math.pi)
    rsoma=L*1e4
    cm_uFcm2=gAIY_scaled[8]

    
    
    soma=h.Section(name="soma")
    soma.L=rsoma
    soma.diam=rsoma
    soma.cm=cm_uFcm2
    soma.Ra=100
    h.psection(sec=soma)
    
    soma.insert('egl19')
    soma.insert('slo1egl19')  
    soma.insert('nca')
    soma.insert('leak')
    soma.insert('slo1iso')
    soma.insert('kqt1')
    soma.insert('shl1')

    
    
    
    for seg in soma:
        
        seg.leak.gbar =  gAIY_scaled[0]
        seg.slo1iso.gbar =  gAIY_scaled[1]
        seg.kqt1.gbar=gAIY_scaled[2]                
        seg.egl19.gbar=gAIY_scaled[3]  
        seg.slo1egl19.gbar = gAIY_scaled[4]        
        seg.nca.gbar =  gAIY_scaled[5]
        seg.shl1.gbar =  gAIY_scaled[6]         
        seg.leak.e=gAIY_scaled[7]
       
    
        seg.eca=60
        seg.ek=-80
        
    stim=h.VClamp(soma(0.5))
    dir(stim)
    
    simdur = 1700
    stim.amp[0]=-60
    stim.amp[2]=-60
    stim.dur[0]=1100
    stim.dur[1]=500
    stim.dur[2]=100
    
    ik_vec = h.Vector()   
    ica_vec=h.Vector()
    inca_vec=h.Vector()
    ileakAIY_vec=h.Vector()
    t_vec = h.Vector()  


   
    
    ik_vec.record(soma(0.5)._ref_ik)
    ica_vec.record(soma(0.5)._ref_ica)
    inca_vec.record(soma(0.5)._ref_i_nca)
    ileakAIY_vec.record(soma(0.5)._ref_i_leak)
 

    t_vec.record(h._ref_t)
    
    ref_ik=[]
    ref_ica=[]
    ref_t=[]
    ref_inca=[]
    ref_ileakAIY=[]
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

               
        #calcium currents
        ref_ica_vec=numpy.zeros_like(ica_vec)
        ica_vec.to_python(ref_ica_vec)
        ref_ica.append(ref_ica_vec)
        
              
        # NCA currents
        ref_inca_vec=numpy.zeros_like(inca_vec)
        inca_vec.to_python(ref_inca_vec)
        ref_inca.append(ref_inca_vec)
    
        # LEAKAGE current
        ref_ileakAIY_vec=numpy.zeros_like(ileakAIY_vec)
        ileakAIY_vec.to_python(ref_ileakAIY_vec)
        ref_ileakAIY.append(ref_ileakAIY_vec)
    
    
    
     
    # total current calculation
    itot=[]
    itot=map(sum, zip(ref_ik,ref_ica,ref_ileakAIY,ref_inca))        
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
    ind=numpy.where(numpy.logical_and(time[0]>=550, time[0]<=599))
    ind_max=numpy.amax(ind)
    ind_min=numpy.amin(ind)
    iv=numpy.mean(itot_normalized[:,ind_min:ind_max],axis=1)
	
	 # CALCULATION OF PEAK CURRENT-VOLTAGE RELATION 
    ind2=numpy.where(numpy.logical_and(time[0]>=100, time[0]<=200))
    ind2_max=numpy.amax(ind2)
    ind2_min=numpy.amin(ind2)
    iv_peak=numpy.amax(itot_normalized[:,ind2_min:ind2_max])
    iv_peak=[]
    
    
    for j in range(ns):
        if j<=6:
            peak=numpy.amin(itot_normalized[j,ind2_min:ind2_max])
        else:
            peak=numpy.amax(itot_normalized[j,ind2_min:ind2_max])
        iv_peak.append(peak)

    return itot_normalized, time, iv_peak, iv

