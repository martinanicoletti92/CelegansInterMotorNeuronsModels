TITLE EGL19
: L-type channels 
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)

}


NEURON {
	SUFFIX egl19
	USEION ca READ eca WRITE ica
    RANGE gbar
	GLOBAL minf, hinf, mtau, htau, megl19,hegl19
	}


PARAMETER{
    v (mV)
    eca (mV)
    celsius (degC)
    gbar=1.55 (S/cm2)
	
    va_egl19=5.6 (mV)
	ka_egl19=7.50 (mV)
	
	shift=10 (mV)
  
    p1hegl19=1.4314 
    p2hegl19=24.8573 (mV)
    p3hegl19=11.9541 (mV)
    p4hegl19=0.1427
    p5hegl19=5.9589
    p6hegl19=-10.5428 (mV)
    p7hegl19=8.0552 (mV)
    p8hegl19=0.6038
  
    pdg1=2.3359 (ms)
	pdg2=2.9324 (ms)
	pdg3=5.2357 (mV)
	pdg4=6.0 (mV)
	pdg5=1.8739 (ms)
	pdg6=1.3930 (mV)
	pdg7=30.0 (mV)
    ctm19=1
     
	pds1=0.4 
	pds2=0.55 
	pds3=81.1179  (ms)
	pds4=-22.9723 (mV)
	pds5=5 (mV)
	pds6=43.0937 (ms)
	pds7=0.9
	pds8=40.4885 (ms)
	pds9=28.7251 (mV)
	pds10=3.7125 (mV)
	pds11=0

	
 
}

ASSIGNED{
		 ica (mA/cm2)
		 minf
		 hinf
		 mtau
		 htau
		 megl19
		 hegl19
	   }
	

STATE {
	m h
}


BREAKPOINT {
	SOLVE states METHOD cnexp
    ica = gbar*m*h*(v-eca)
	megl19=m
    hegl19=h	
}


INITIAL {
	rates(v)
	m=minf
	h=hinf
	
}

DERIVATIVE states {    
     rates(v) 
        m' = (minf - m)/mtau
		h'=(hinf-h)/htau
		
	}



FUNCTION act(v (mV)) { 
		
        act=1/(1+exp(-(v-va_egl19+shift)/ka_egl19))
		
		}
FUNCTION inact(m (mV)){
		
        inact = ((p1hegl19/(1+exp(-(v-p2hegl19+shift)/p3hegl19))+p4hegl19)*(p5hegl19/(1+exp((v-p6hegl19+shift)/p7hegl19))+p8hegl19))
        
		}
FUNCTION tact(m (mV)){
		
		tact= (pdg1+(pdg2*exp(-(v-pdg3+shift)^2/(pdg4)^2))+(pdg5*exp(-(v-pdg6+shift)^2/(pdg7)^2)))*ctm19
		
		}
		
FUNCTION tinact(v(mV)){
		
		tinact= pds1*(((pds2*pds3)/(1+exp((v-pds4+shift)/pds5)))+pds6+((pds7*pds8)/(1+exp((v-pds9+shift)/pds10)))+pds11)
		
		}

PROCEDURE rates(v (mV)){
       minf=act(v)
	   hinf=inact(v)
	   mtau=tact(v)
       htau=tinact(v)
}

















