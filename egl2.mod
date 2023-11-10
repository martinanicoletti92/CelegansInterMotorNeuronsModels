TITLE EGL2
: egl2 currents model
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
	(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX egl2
	USEION k READ ek WRITE ik
    RANGE gbar
	
}


PARAMETER{
   v (mV)
   ek (mV)
   celsius		(degC)
   gbar=0.85 (S/cm2)
   va_egl2=-6.8594 (mV)
   ka_egl2=14.9131 (mV)
   stmegl2=0  
   cegl2=0.5
   
  p1tmegl2=16.7800
  p2tmegl2=-122.5682
  p3tmegl2=13.7976
  p4tmegl2=8.0969
  fegl2=1
}

ASSIGNED{
		 ik (mA/cm2)
		
    }
	

STATE {
	m 
}


BREAKPOINT {
	SOLVE states METHOD cnexp

	ik = gbar*m*(v-ek)

}


INITIAL {
		m=minf(v)
	
}

DERIVATIVE states {     
        
        m' = (minf(v) - m)/mtau(v)
		
	}



FUNCTION minf(v (mV)) {
		
        minf=1/(1+exp(-(v-va_egl2+stmegl2)/(ka_egl2*fegl2)))
		
		}
		
FUNCTION mtau(v (mV)){
		
        mtau=((p1tmegl2/(1+exp((v-p2tmegl2+stmegl2)/p3tmegl2)))+p4tmegl2)*cegl2
		
		}




