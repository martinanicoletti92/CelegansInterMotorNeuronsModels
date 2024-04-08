: TITLE SLO1_ISO
: isolated SLO1 channels 
: "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
: M. Nicoletti et al. PloS ONE, 19(3): e0298105.
: https://doi.org/10.1371/journal.pone.0298105


NEURON {
       SUFFIX slo1iso
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gbar, c1


}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
       v (mV)
	   cai (mM)
	   fondo=0.05 (uM)
	   ek=-80 (mV)
	   eca=60 (mV)
	  
	   celsius	(degC)
	   gbar=.11 (S/cm2)
	   wom=3.152961  (/ms)  
	   wyx=0.012643 (/mV)
	   kyx=34.338784 (uM)
	   nyx=0.000100   (1)
	   wop=0.156217 (/ms)
	   wxy=-0.027527 (/mV)
	   kxy=55.726186 (/ms)
	   nxy=1.299198 (1) 
	   c1=1 (1)
   
	
}

ASSIGNED {
         minf
         mtau          (ms)
		 s0 (mV)
		 v0 (mV)
         ik            (mA/cm2)
}

STATE {
      m   FROM 0 TO 1

}

BREAKPOINT {
           SOLVE states METHOD cnexp
           ik = gbar * m * (v - ek)
}

DERIVATIVE states {
        rates(v, cai)

        m' = (minf - m) / mtau

}

PROCEDURE rates(v (mV), ca (mM)) {


		UNITSOFF
		  s0=1/(wyx-wxy)
		  v0=s0*(log(wom/wop)+log(1+pow(kxy/(ca*1e3),nxy))-log(1+pow((ca*1e3)/kyx,nyx)))
          minf =1/(1+exp(-(v-v0)/s0))
          mtau = ((exp(wxy*v)/wop)*(1+pow(kxy/(ca*1e3),nxy))*(1/(1+exp(-(v-v0)/s0))))*c1
	    UNITSON
}

INITIAL {
        rates(v,cai)
        m = minf
}
