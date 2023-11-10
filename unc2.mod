TITLE UNC2
: P/Q-type channels 
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738) 



UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)

}


NEURON {
	SUFFIX unc2
	USEION ca READ eca WRITE ica
    RANGE gbar
	GLOBAL minf, hinf, mtau, htau, munc2,hunc2
	}


PARAMETER{
    v (mV)
    eca (mV)
    celsius (degC)
    gbar=1 (S/cm2)
    va_unc2=-12.17 (mV)
	ka_unc2=3.97 (mV)
	vi_unc2=-52.47 (mV)
	ki_unc2=5.6 (mV)
	stm2=25 (mV)
	sth2=25 (mV)
   
    p1tmunc2=1.4969 (ms)
	p2tmunc2=-8.1761 (mV)
	p3tmunc2=9.0753 (mV)
	p4tmunc2=15.3456 (mV)
	p5tmunc2=0.1029 (ms)
  
    p1thunc2=83.8037 (ms)
	p2thunc2=52.8997 (mV)
	p3thunc2=3.4557 (mV)
	p4thunc2=72.0995 (ms)
	p5thunc2=23.9009 (mV)
	p6thunc2=3.5903 (mV)
	
	fp3=1
	fp4=1
	fp5=1
	
	shifthunc2=30
	shiftmunc2=30
	consthunc2=1.7
	constmunc2=3
	func2=1
	f2unc2=1
}

ASSIGNED{
		 ica (mA/cm2)
		 minf
		 hinf
		 mtau
		 htau
		 munc2
		 hunc2
	   }

STATE {
	m h
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gbar*m*h*(v-eca)
	munc2=m
    hunc2=h	
}


INITIAL {
	rates(v)
	m=minf
	h=hinf
	
}

DERIVATIVE states {    
     rates(v) 
        m' = (minf - m)/mtau
		h'=(hinf-h)/htau
		
	}


FUNCTION act(v(mV)){

     act=1/(1+exp(-(v-va_unc2+stm2)/ka_unc2))

}

FUNCTION inact(v(mV)){

     inact= 1/(1+exp((v-vi_unc2+sth2)/ki_unc2))

}


FUNCTION tact(v(mV)){
      
      tact=(p1tmunc2/(exp(-(v-p2tmunc2+shiftmunc2)/(p3tmunc2*fp3))+exp((v-p2tmunc2+shiftmunc2)/(p4tmunc2*fp4)))+p5tmunc2)*constmunc2

	  }

FUNCTION tinact(v(mV)){
   
      tinact=(p1thunc2/(1+exp((v-p2thunc2+shifthunc2)/(p3thunc2*fp5)))+p4thunc2/(1+exp(-(v-p5thunc2+shifthunc2)/(p6thunc2*fp5))))*consthunc2
   
	  }
	  
	  
PROCEDURE rates(v (mV)){
       minf=act(v)
	   hinf=inact(v)
	   mtau=tact(v)
       htau=tinact(v)
}
	  
	  

