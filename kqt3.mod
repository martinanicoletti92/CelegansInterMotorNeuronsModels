TITLE KQT3
: kqt3 channels
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}


NEURON {
	SUFFIX kqt3
	USEION k READ ek WRITE ik
    RANGE gbar
	}


PARAMETER{
   v (mV)
   ek (mV)
   celsius (degC)
   gbar=0.55 (S/cm2)
   
   va_kqt3=-12.6726 (mV)
   ka_kqt3=15.8008 (mV)
  
     
   p1tmskqt3=5503 (ms)
   p2tmskqt3=5345.4 (ms)
   p3tmskqt3=-0.02827 (1/mV)
   p4tmskqt3=-23.9 (mV)
   p5tmskqt3=4590.6 (ms)
   p6tmskqt3=-0.0357 (1/mV)
   p7tmskqt3=14.15   (mV)
  
   p1tmfkqt3=395.3 (ms)
   p2tmfkqt3=38.1 (mV)
   p3tmfkqt3=33.59 (mV)
   constkqt3=10 
   
   w1=0.49
   w2=0.51
   w3=1.084
   w4=28.78
   
   tw1=5.44 (ms)
   tw2=29.2 (ms)
   tw3=48.09 (mV)
   tw4=48.83 (mV)
    
   sq1=0.34
   sq2=0.66
   sq3=45.3
   sq4=12.3
   tsq1=5000
   ckqt3=0.1
  
}

ASSIGNED{
		 ik (mA/cm2)
		 curr (mA/cm2)
		 g (S/cm2)
	
}	

STATE {
	mf ms s w
}


BREAKPOINT {
	SOLVE states METHOD cnexp

	ik = gbar*((0.3*mf)+(0.7*ms))*s*w*(v-ek)

}


INITIAL {
	
	mf=minf(v)
	ms=minf(v)
    w=winf(v)
	s=sinf(v)
}

DERIVATIVE states {     
        mf' = (minf(v) - mf)/mtauf(v)
		ms'=(minf(v)-ms)/mtaus(v)
		s'=(sinf(v)-s)/500
		w'=(winf(v)-w)/tw(v)

	}

FUNCTION minf(v (mV)) { 
        minf=1/(1+exp(-(v-va_kqt3+constkqt3)/ka_kqt3))
	}


FUNCTION mtauf(v (mV)){

        mtauf=(p1tmfkqt3/(1+((v+p2tmfkqt3)/p3tmfkqt3)^2))*ckqt3

		}
FUNCTION mtaus(v (mV)){

		mtaus=(p1tmskqt3-p2tmskqt3/(1+10^(p3tmskqt3*(p4tmskqt3-v)))-p5tmskqt3/(1+10^(p6tmskqt3*(v+p7tmskqt3))))*ckqt3
	}

FUNCTION winf(v (mV)){
	winf=w1+(w2/(1+exp((v+w3)/w4)))
}

FUNCTION sinf(v (mV)){
	sinf=sq1+(sq2/(1+exp((v+sq3)/sq4)))
}

FUNCTION tw(v (mV)){
	tw=(tw1+(tw2/(1+pow((v+tw3)/tw4,2))))*ckqt3
}
		










