TITLE SHL1
: Shl1 current cultured myocites
: Fawcett et. al 2006 https://doi.org/10.1074/jbc.M605814200
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
	(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX shl1
	USEION k READ ek WRITE ik
    RANGE gbar
	}


PARAMETER{
   v (mV)
   ek (mV)
   celsius	(degC)
   gbar=2.9 (S/cm2)
   vashal=10.26 (mV)
   kashal=16.250 (mV)
   kishal=8.3 (mV)
   vishal=-40 (mV)
   shalsfhit=10 (mV)

   
   ptmshal1=13.8  (ms)
   ptmshal2=-40 (mV)
   ptmshal3=12.9213 (mV)
   ptmshal4=-40 (mV)
   ptmshal5=6.4876 (mV)
   ptmshal6=1.8849 (ms)
  
   pthfshal1=539.1584  (ms)
   pthfshal2=-60 (mV)
   pthfshal3=4.9199   (mV)
   pthfshal4=27.2811  (ms)

   pthsshal1=8422     (ms)
   pthsshal2=-60 (mV)
   pthsshal3=6.3785    (mV)
   pthsshal4=118.8983  (ms)
   
   a=0.8

  
}

ASSIGNED{
		 ik (mA/cm2)
		
    }
	

STATE {
	m hs hf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*m*((a*hf)+((1-a)*hs))*(v-ek)

}


INITIAL {
	m=minf(v)
	hs=hinf(v)
	hf=hinf(v)
}

DERIVATIVE states {     
        m' = (minf(v) - m)/mtau(v)
		hs'=(hinf(v)-hs)/htaus(v)
		hf'=(hinf(v)-hf)/htauf(v)
}


FUNCTION minf(v (mV)) { 
        UNITSOFF	
        minf=1/(1+exp(-(v-vashal+shalsfhit)/kashal)) 
	    UNITSON
		}
		
FUNCTION hinf(v (mV)){
		UNITSOFF
        hinf = 1/(1+exp((v-vishal+shalsfhit)/kishal))
		UNITSON
		}
		
FUNCTION mtau(v(mV)){
	    UNITSOFF
        mtau=(ptmshal1/(exp(-(v-ptmshal2)/ptmshal3)+exp((v-ptmshal4)/ptmshal5))+ptmshal6)/2
		UNITSON
		}
		
FUNCTION htauf(v(mV)){
	    UNITSOFF
		htauf=(pthfshal1/(1+exp((v-pthfshal2)/pthfshal3))+pthfshal4)/3
		UNITSON
		}
		
FUNCTION htaus(v(mV)){
		
		htaus=(pthsshal1/(1+exp((v-pthsshal2)/pthsshal3))+pthsshal4)
		
}
