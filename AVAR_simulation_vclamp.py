# "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
# M. Nicoletti et al. PloS ONE, 19(3): e0298105.
# https://doi.org/10.1371/journal.pone.0298105

def AVA_simulation_vc(gAVA_scaled,vstart,vstop,ns):
    
    # voltage clamp per Minimizzazioni con canali calcio fissati e ottimizzazione sul potassio
    #i valori delle conduttanze del calcio vengono dal set xpp DEFINITIVO.ode
    #fissati anche i valori di irk, kcnl per fittare le code e di slo1 slo2
    #libere anche le conduttanze di perdita
 
    from neuron import h,gui
    import numpy
    import math
   
    
    surf=1121.79e-8 # surface in cm^2 form neuromorpho AVAL
    vol=129.08e-12 # total volume
    Fc=96485 #Faraday constant C/mol
    Z=2 #Valence
    alphan=surf/(Z*Fc*vol)
    L=math.sqrt(surf/math.pi)
    rsoma=L*1e4

    cm_uFcm2=gAVA_scaled[6]
    cm_pF=cm_uFcm2*1e6/surf

    #fixed conductances
   
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
    soma.insert('unc103')

    for seg in soma:
        

        seg.egl19.gbar=gAVA_scaled[0]
        seg.leak.gbar=gAVA_scaled[1]
        seg.irk.gbar=gAVA_scaled[2]
        seg.nca.gbar=gAVA_scaled[3]
        seg.unc103.gbar=gAVA_scaled[4]      
        seg.leak.e=gAVA_scaled[5]
        
        seg.eca=60
        seg.ek=-80
        
    stim=h.VClamp(soma(0.5))
    dir(stim)
    
    simdur = 1500
    stim.amp[0]=-30
    stim.amp[2]=-30
    stim.dur[0]=1007.8
    stim.dur[1]=250
    stim.dur[2]=242.2
    
    ik_vec = h.Vector()   
    ica_vec=h.Vector()
    inca_vec=h.Vector()
    
    ileakAVA_vec=h.Vector()
    t_vec = h.Vector()  
    inca_vec.record(soma(0.5)._ref_i_nca)
    ik_vec.record(soma(0.5)._ref_ik)
    ica_vec.record(soma(0.5)._ref_ica)    
    ileakAVA_vec.record(soma(0.5)._ref_i_leak)


    t_vec.record(h._ref_t)
    
    ref_ik=[]
    ref_ica=[]
    ref_inca=[]
    ref_t=[]

    ref_ileakAVA=[]

        
    for i in numpy.linspace(start=vstart, stop=vstop, num=ns):
        
        stim.amp[1]=i
        h.tstop=simdur
        h.dt=0.01
        h.finitialize(-30)
        h.run()
        
        #time
        ref_t_vec=numpy.zeros_like(t_vec)
        t_vec.to_python(ref_t_vec)
        ref_t.append(ref_t_vec)
        
        # potassium current
        ref_ik_vec=numpy.zeros_like(ik_vec)
        ik_vec.to_python(ref_ik_vec)
        ref_ik.append(ref_ik_vec)
        
        # NCA currents
        ref_inca_vec=numpy.zeros_like(inca_vec)
        inca_vec.to_python(ref_inca_vec)
        ref_inca.append(ref_inca_vec)
               
        #calcium currents
        ref_ica_vec=numpy.zeros_like(ica_vec)
        ica_vec.to_python(ref_ica_vec)
        ref_ica.append(ref_ica_vec)


        # LEAKAGE current
        ref_ileakAVA_vec=numpy.zeros_like(ileakAVA_vec)
        ileakAVA_vec.to_python(ref_ileakAVA_vec)
        ref_ileakAVA.append(ref_ileakAVA_vec)
    
    
    
     
    # total current calculation
    itot=[]
    itot=map(sum, zip(ref_ik,ref_ica,ref_inca,ref_ileakAVA))         
    current=numpy.array(list(itot))
    inorm=current*1e9*surf #total current in pA
    #time array
    time1=numpy.array(ref_t)
    resc_ind=numpy.where(time1[1,:]>=1000)
    resc_min=numpy.amin(resc_ind)
    resc_max=numpy.amax(resc_ind)
    itot_normalized=inorm[:,resc_min:resc_max]
    time=time1[:,resc_min:resc_max]-1000
    
    
    ## SS I-V CURVE
    ind=numpy.where(numpy.logical_and(time[0]>=247, time[0]<=257))
    ind_max=numpy.amax(ind)
    ind_min=numpy.amin(ind)
    iv=numpy.mean(itot_normalized[:,ind_min:ind_max],axis=1)
	
	 # PEAKS I-V CURVE
    ind2=numpy.where(numpy.logical_and(time[0]>=7.8, time[0]<=17.8))
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



