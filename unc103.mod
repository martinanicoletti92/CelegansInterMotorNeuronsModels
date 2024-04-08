TITLE UNC-103
: "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
: M. Nicoletti et al. PloS ONE, 19(3): e0298105.
: https://doi.org/10.1371/journal.pone.0298105


UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX unc103
	USEION k READ ek WRITE ik
    RANGE gbar,g,curr
	}


PARAMETER{
   v (mV)
   ek (mV)
   celsius	(degC)
   gbar=2.9 (S/cm2)
   va=-15.1 (mV)
   ka=7.85 (mV)
   vi=-48 (mV)
   ki=28 (mV)
   th1=8.1559 (ms)
   th2=-25.2890 (mV)
   th3=29.5074 (mV)
   th4=0.2300 (ms)
   tm1=87.4088 (ms)
   tm2=-28.3339 (mV)
   tm3=13.0998 (mV)
   tm4=0.2562 (ms)
}

ASSIGNED{
		 ik (mA/cm2)
		 curr (mA/cm2)
		 g (S/cm2)
		
    }
	

STATE {
	m h
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g=gbar*m*h
	curr=gbar*m*h*(v-ek)
	ik = gbar*m*h*(v-ek)

}


INITIAL {
	m=minf(v)
	h=hinf(v)
}

DERIVATIVE states {     
        m' = (minf(v) - m)/mtau(v)
		h'=(hinf(v)-h)/htau(v)
}


FUNCTION minf(v (mV)) { 
        UNITSOFF	
        minf=1/(1+exp(-(v-va)/ka)) 
	    UNITSON
		}
		
FUNCTION hinf(v (mV)){
		UNITSOFF
        hinf = 1/(1+exp((v-vi)/ki))
		UNITSON
		}
		
FUNCTION mtau(v(mV)){
	    UNITSOFF
        mtau=((tm1/(1+exp((v-tm2)/tm3)))+tm4)*((tm1/(1+exp(-(v-tm2)/tm3)))+tm4)
		UNITSON
		}
		
FUNCTION htau(v(mV)){
	    UNITSOFF
		htau=((th1/(1+exp((v-th2)/th3)))+th4)*((th1/(1+exp(-(v-th2)/th3)))+th4)
		UNITSON
		}
