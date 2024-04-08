TITLE KQT1
: "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
: M. Nicoletti et al. PloS ONE, 19(3): e0298105.
: https://doi.org/10.1371/journal.pone.0298105


UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX kqt1
	USEION k READ ek WRITE ik
    RANGE gbar,g,curr
	}


PARAMETER{
   v (mV)
   ek (mV)
   celsius	(degC)
   gbar=2.9 (S/cm2)
   va=-17.6053 (mV)
   ka=9.5843 (mV)
   p1tmkqt1=10 (ms)
   p2tmkqt1=895.9 (ms)
   p3tmkqt1=-18.01 (mV)
   p4tmkqt1=31.04 (mV)
   s1=0.41 
   s2=-86.84 (mV)
   s3=15.05 (mV)
   s4=0.59 
   s5=70.13 (mV)
   s6=13.37 (mV)
   p1tskqt1=1077 (ms)
   p2tskqt1=185845 (ms)
   p3tskqt1=39.44 (mV)
   p4tskqt1=7.34 (mV)
}





ASSIGNED{
		 ik (mA/cm2)
		 curr (mA/cm2)
		 g (S/cm2)
		
    }
	

STATE {
	m s
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g=gbar*m
	curr=gbar*m*s*(v-ek)
	ik = gbar*m*s*(v-ek)

}

INITIAL {
	m=minf(v)
	s=sinf(v)
}

DERIVATIVE states {     
        m' = (minf(v) - m)/mtau(v)
	    s' = (sinf(v) - s)/stau(v)
}


FUNCTION minf(v (mV)) { 
        UNITSOFF	
        minf=1/(1+exp(-(v-va)/ka)) 
	    UNITSON
		}
		

FUNCTION sinf(v (mV)) { 
        UNITSOFF	
        sinf=(s1/(1+exp((v-s2)/s3)))+(s4/(1+exp((v-s5)/s6)))
	    UNITSON
		}		
		
		
FUNCTION mtau(v(mV)){
	    UNITSOFF
        mtau=(p2tmkqt1/(1+exp((p3tmkqt1-v)/p4tmkqt1)))+p1tmkqt1
		UNITSON
		}
		
FUNCTION stau(v(mV)){
	    UNITSOFF
        stau=(p1tskqt1+(p2tskqt1/(1+((v-p3tskqt1)/p4tskqt1)^2)))
		UNITSON
		}
