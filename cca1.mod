TITLE CCA1
: T-type, CCA-1, channels 
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)


UNITS {
	(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)

}


NEURON {
	SUFFIX cca1
	USEION ca READ eca WRITE ica
    RANGE gbar
	}


PARAMETER{
    v (mV)
    eca (mV)
    celsius (degC)
    gbar=0.7 (S/cm2)
    va_cca1=-42.65
	ka_cca1=1.7
	sscca1=15
	stmcca1=30
	sthcca1=15
	sshcca1=15
	constmcca1=0.5
	consthcca1=0.08
	fcca=1.4
	f2cca1=1.15
	f3ca=1.7
	f4ca=1.1
	vi_cca1=-58
	ki_cca1=7
   
    p1tmcca1=40
	p2tmcca1=-62.5393
	p3tmcca1=-12.4758
	p4tmcca1=0.6947
  
    p1thcca1=280
	p2thcca1=-60.7312
	p3thcca1=8.5224
	p4thcca1=19.7456
	
	
}

ASSIGNED{
		 ica (mA/cm2)
		}
	

STATE {
	m h
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gbar*m*m*h*(v-eca)

}


INITIAL {
	m=minf(v)
	h=hinf(v)
	
}

DERIVATIVE states {     
        m' = (minf(v) - m)/mtau(v)
		h'=(hinf(v)-h)/htau(v)
	}


FUNCTION minf(v (mV)){
 
       minf=1/(1+exp(-(v-va_cca1+sscca1)/(ka_cca1*fcca)))

}

FUNCTION hinf(m (mV)){
	
	hinf =1/(1+exp((v-vi_cca1+sshcca1)/(ki_cca1*f2cca1)))
	
}

FUNCTION mtau(v (mV)){
	
	mtau=((p1tmcca1/(1+exp(-(v-p2tmcca1+stmcca1)/(p3tmcca1*f3ca))))+p4tmcca1)*constmcca1
	
}

FUNCTION htau(v (mV)) { 
		
        htau=((p1thcca1/(1+exp((v-p2thcca1+sthcca1)/(p3thcca1*f4ca))))+p4thcca1)*consthcca1
		
		}
