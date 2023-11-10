TITLE KVS1
: kvs1
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)

}


NEURON {
	SUFFIX kvs1
	USEION k READ ek WRITE ik
    RANGE gbar
	
}


PARAMETER{
   v (mV)
   ek (mV)
   celsius		(degC)
   gbar=0.8 (S/cm2)
   va_kvs1=57.1 (mV)
   ka_kvs1=25   (mV)
   vi_kvs1=47.3  (mV)
   ki_kvs1=11.1   (mV)
   
   p1tmkvs1=30.0000 (ms)
   p2tmkvs1=18.1232 (mV)
   p3tmkvs1=-20.0000 (mV)
   p4tmkvs1=1.000    (ms)
  
   p1thkvs1=88.4715 (ms)
   p2thkvs1=50.00   (mV)
   p3thkvs1=-15      (mV)
   p4thkvs1=53.4060  (ms)
   skvs1=30           (mV)
   cthkvs1=0.1       
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
        minf=1/(1+exp(-(v-va_kvs1+skvs1)/ka_kvs1))
		UNITSON
		}
		
FUNCTION hinf(v (mV)){
		UNITSOFF
        hinf = 1/(1+exp((v-vi_kvs1+skvs1)/(ki_kvs1)))
		UNITSON
		}
		
FUNCTION mtau(v (mV)){
		UNITSOFF
        mtau=(p1tmkvs1/(1+exp(-(v-p2tmkvs1)/p3tmkvs1))+p4tmkvs1)*cthkvs1
		UNITSON
		}

FUNCTION htau(v (mV)){
		UNITSOFF
		htau=(p1thkvs1/(1+exp(-(v-p2thkvs1)/p3thkvs1))+p4thkvs1)*cthkvs1
		UNITSON
		}












