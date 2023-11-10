TITLE IRK
: IRK current model
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
	(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)

}


NEURON {
	SUFFIX irk
	USEION k READ ek WRITE ik
    RANGE gbar
	}


PARAMETER{
   v (mV)
   ek (mV)
   celsius		(degC)
   gbar=.65 (nS/cm2)
   va_kir=-52 (mV)
   ka_kir=13  (mV)
   
   p1tmkir=17.0752
   p2tmkir=-17.8258
   p3tmkir=20.3154
   p4tmkir=-43.4414
   p5tmkir=11.1691
   p6tmkir=3.8329
  
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
      
        minf=1/(1+exp((v-va_kir+30)/ka_kir))
	  
	  }
	  
FUNCTION mtau(v(mV)){
		
        mtau=p1tmkir/(exp(-(v-p2tmkir)/p3tmkir)+exp((v-p4tmkir)/p5tmkir))+p6tmkir
		
		}
		












