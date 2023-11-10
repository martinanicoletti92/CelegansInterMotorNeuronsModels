TITLE EGL36
: egl36 currents model 
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
	(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)

}


NEURON {
	SUFFIX egl36
	USEION k READ ek WRITE ik
    RANGE gbar
	
}


PARAMETER{
   v (mV)
   ek (mV)
   celsius		(degC)
   gbar=0.8 (S/cm2)
   va_egl36=63 (mV)
   ka_egl36=28.5 (mV)
   
   t1_egl36=355 (ms)
   t2_egl36=63  (ms)
   t3_egl36=13  (ms)
   
   a1=0.31     (1)
   a2=0.36     (1)
   a3=0.39     (1)
 
}

ASSIGNED{
		 ik (mA/cm2)
		 
    }
	

STATE {
	m1 m2 m3 
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*((a1*m1)+(a2*m2)+(a3*m3))*(v-ek)

}


INITIAL {
	m1=minf(v)
	m2=minf(v)
	m3=minf(v)
	
}

DERIVATIVE states {     
        m1' = (minf(v) - m1)/t1_egl36
		m2' = (minf(v) - m2)/t2_egl36
		m3' = (minf(v) - m3)/t3_egl36
	}


      
FUNCTION minf(v (mV)) {
		
        minf=1/(1+exp(-(v-va_egl36)/ka_egl36))
		
		}



