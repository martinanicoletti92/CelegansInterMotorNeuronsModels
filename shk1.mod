TITLE SHK1
: Shk1 current cultured myocites
: Fawcett 2006 inactivation and tau (https://doi.org/10.1074/jbc.M605814200)
: activation from Liu et al 2018 (DOI: 10.1016/j.cell.2018.08.018)
: From Nicoletti et al 2023



UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX shk1
	USEION k READ ek WRITE ik
    RANGE gbar
	}


PARAMETER{
   v (mV)
   ek (mV)
   celsius		(degC)
   gbar=0.1 (S/cm2)
   vashak=2
   kashak=10
   kishak=5.8
   vishak=-6.95
   
   ptmshak1=26.571450568169027 (ms)
   ptmshak2=-33.741611800716130 (mV)
   ptmshak3=15.757936311607475  (mV)
   ptmshak4=15.364937728953288 (mV)
   ptmshak5=1.990037272604829 (ms)
   shiftV05=0
  
  pthshak=1400 (ms)
  
  
}

ASSIGNED{
		 ik (mA/cm2)
			 
    }
	

STATE {
	m h
}


BREAKPOINT {
	SOLVE states METHOD cnexp
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
        minf=1/(1+exp(-(v-vashak)/(kashak)))
		UNITSON
		}
		
FUNCTION hinf(v(mV)){
		UNITSOFF
        hinf = 1/(1+exp((v-vishak)/(kishak)))
		UNITSON
		}

FUNCTION mtau (v(mV)){	
		UNITSOFF
		mtau=ptmshak1/(exp(-(v-(ptmshak2+shiftV05))/ptmshak4)+exp((v-(ptmshak2+shiftV05))/ptmshak3))+ptmshak5
		UNITSON
		}
		
FUNCTION htau(v(mV)){
	    UNITSOFF
		htau=pthshak
		UNITSON
}












