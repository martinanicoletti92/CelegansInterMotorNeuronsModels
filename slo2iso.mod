: TITLE SLO2ISO
: isolated SLO2 channels 
: "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
: M. Nicoletti et al. PloS ONE, 19(3): e0298105.
: https://doi.org/10.1371/journal.pone.0298105


NEURON {
       SUFFIX slo2iso
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gbar, ik,c2
       GLOBAL minf, mtau
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
	   v (mV)
	   cai (mM)
	   ek=-80 (mV)
	   eca=60 (mV)

	  
	   celsius		(degC)
	   gbar=.1 (S/cm2)
	   wom1=0.896395  (/ms)  
	   wyx1=0.019405  (/mV)
	   kyx1=3294.553404 (uM)
	   nyx1=0.000010   (1)
	   wop1=0.026719 (/ms)
	   wxy1=-0.024123 (/mV)
	   kxy1=93.449423 (/ms)
	   nxy1=1.835067 (1)
	   c2=1 (1)
	   
  
	
}

ASSIGNED {
         minf
         mtau          (ms)
		 s0 (mV)
		 v0 (mV)
         ik  (mA/cm2)
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
		  s0=1/(wyx1-wxy1)
		  v0=s0*(log(wom1/wop1)+log(1+pow(kxy1/(ca*1e3),nxy1))-log(1+pow((ca*1e3)/kyx1,nyx1)))
          minf =1/(1+exp(-(v-v0)/s0))
          mtau = ((exp(wxy1*v)/wop1)*(1+pow(kxy1/(ca*1e3),nxy1))*(1/(1+exp(-(v-v0)/s0))))*c2
		UNITSON
}

INITIAL {
        rates(v,cai)
        m = minf
}
