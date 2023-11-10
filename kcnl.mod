TITLE sk
: kcnl channels
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)



UNITS {
     (mA) = (milliamp)
	 (S) = (siemens)
	 (mV) = (millivolt)
	 (molar)=(1/liter)
	 (uM) =	(micromolar)
	 (mM)= (millimolar)
	 (M)= (molar)
	
}

NEURON {
	SUFFIX kcnl
	USEION k READ ek WRITE ik
	USEION ca READ eca
	RANGE  gbar, g, curr
	EXTERNAL calcium_caintra1
 
}


PARAMETER{
   v (mV)
   eca (mV)
   ek (mV)
   celsius		(degC)
   ca_eq=0.05e-6 (M)
   gbar=0.06 (S/cm2)
   km=0.33e-6 (M)
   t_sk=6.3 (ms)
      
}

ASSIGNED{
		 ik (mA/cm2)
		 curr (mA/cm2)
		 cai (M)
		 minf
         g (S/cm2)	
         		 
		
}	

STATE {
	m
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	g=m
	curr=gbar*m*(v-ek)
	ik = gbar*m*(v-ek)

}


INITIAL {
	
	m=act(ca_eq)
	
}

DERIVATIVE state{  
        
        m' = (act(calcium_caintra1) - m)/t_sk
	
	}


FUNCTION act(ca (M)) { 
	UNITSOFF
        act=(ca)/((ca)+km)
	UNITSON
	}
	


