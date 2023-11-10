TITLE EXP2
: From Nicoletti et al 2023

UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX exp2
	USEION k READ ek WRITE ik
    RANGE gbar
	}


PARAMETER{
   v (mV)
   ek (mV)
   gbar (S/cm2)
   celsius		(degC)
   va_exp2=-17 (mV)
   ka_exp2=6.5 (mV)
   vi_exp2=-50 (mV)
   ki_exp2=10 (mV)
   
  ptm1_exp2=40
  ptm2_exp2=-20.5284
  ptm3_exp2=4.7071
  ptm4_exp2=2.4025
  pth1_exp2=1.2790
  pth2_exp2=-89.0175
  pth3_exp2=49.3492
  pth4_exp2=0.8739
 
  
}

ASSIGNED{
		 ik (mA/cm2)
		 g(S/cm2)
		 curr(mA/cm2)
			 
    }
	

STATE {
	m h
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*m*h*h*(v-ek)

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
        minf=1/(1+exp(-(v-va_exp2)/ka_exp2))
		UNITSON
		}
		
FUNCTION hinf(v(mV)){
		UNITSOFF
        hinf = 1/(1+exp((v-vi_exp2)/ki_exp2))
		UNITSON
		}

FUNCTION mtau (v(mV)){	
		UNITSOFF
		mtau=((ptm1_exp2/(1+exp((v-ptm2_exp2)/ptm3_exp2)))+ptm4_exp2)*((ptm1_exp2/(1+exp(-(v-ptm2_exp2)/ptm3_exp2)))+ptm4_exp2)
		UNITSON
		}
		
FUNCTION htau(v(mV)){
	UNITSOFF
		htau=((pth1_exp2/(1+exp((v-pth2_exp2)/pth3_exp2)))+pth4_exp2)*((pth1_exp2/(1+exp(-(v-pth2_exp2)/pth3_exp2)))+pth4_exp2)
		UNITSON
}