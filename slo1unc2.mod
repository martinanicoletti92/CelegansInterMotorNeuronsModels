TITLE slo1unc2
: slo1 channels coupled with unc2 calcium channels (1:1 stoichiometry)
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)

UNITS {
	(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	 (molar)=(1/liter)
	 (uM) =	(micromolar)
	 FARADAY = (faraday) (coulombs)
	 
}

NEURON {
	SUFFIX slo1unc2
	USEION k READ ek WRITE ik
	USEION ca READ eca
	RANGE  gbar
	EXTERNAL  munc2_unc2, hunc2_unc2
 }


PARAMETER{
   v (mV)
   cai (uM)
   bkg=0.05 (uM)
   ek (mV)
   eca (mV)
   munc2_unc2
   hunc2_unc2
  
   celsius		(degC)
   gbar=.11 (S/cm2)
   wom=3.152961  (/ms)  
   wyx=0.012643 (/mV)
   kyx=34.338784 (uM)
   nyx=0.000100   (1)
   wop=0.156217 (/ms)
   wxy=-0.027527 (/mV)
   kxy=55.726186 (/ms)
   nxy=1.299198 (1)
   
   r=13e-9 (nm)
   d=250e-12 (um2/s)
   kb=500e6 (/M-s)
   b=30e-6 (M)
   gsc=40e-12 (S)
   pi=3.14
   
   
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
	
	shifthunc2=30
	shiftmunc2=30
	consthunc2=1.7
	constmunc2=3
	func2=1
	f2unc2=1
	fp3=1
	fp4=1
	fp5=1
    }

ASSIGNED{
		 ik (mA/cm2)
		  g (S/cm2)
		  curr (mA/cm2)
		 minf
		 tslo1
		 alpha1
         beta1	
         ts(ms)	
         v1
         v2	
         ta (ms)
         ti (ms)

}	


STATE {
	m 
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*m*hunc2_unc2*(v-ek)
}


INITIAL {
	rates(calcium(v), v)
	m=minf
}

DERIVATIVE states { 
       rates(calcium(v), v)   
        m' = (minf - m)/tslo1
	
	}


PROCEDURE rates(calcium(v),v (mV)){ 
    alpha1=minfUNC2(v)/tmUNC2(v)
	beta1=(1/tmUNC2(v))-alpha1
 	minf=(munc2_unc2*kop(calcium(v),v)*(alpha1+beta1+kcm(v)))/((kop(calcium(v),v)+kom(calcium(v),v))*(kcm(v)+alpha1)+(beta1*kcm(v)))
	tslo1=((alpha1+beta1+kcm(v))/((kop(calcium(v),v)+kom(calcium(v),v))*(kcm(v)+alpha1)+(beta1*kcm(v))))
 }

		
FUNCTION kcm(v (mV)){
		kcm=wom*exp(-wyx*v)*(1/(1+((bkg/kyx)^nyx)))
}
	
FUNCTION kom(calcium(v),v (mV)){
		kom=wom*exp(-wyx*v)*(1/(1+pow(calcium(v)/kyx,nyx)))
}	

FUNCTION kop(calcium(v),v (mV)){
		kop=wop*exp(-wxy*v)*(1/(1+pow(kxy/calcium(v),nxy)))
}

FUNCTION calcium(v (mV)){
	calcium=(((fabs(gsc*(v-eca)*1e-3)/(8*pi*r*d*FARADAY))*exp(-r/sqrt(d/(kb*b))))*1e6*1e-3)+bkg
}


FUNCTION minfUNC2(v(mV)){

UNITSOFF
     minfUNC2=1/(1+exp(-(v-va_unc2+stm2)/(ka_unc2*func2)))
  UNITSON

}

FUNCTION hinfUNC2(v(mV)){
  UNITSOFF
     hinfUNC2= 1/(1+exp((v-vi_unc2+sth2)/(ki_unc2*f2unc2)))
  UNITSON
}

FUNCTION tmUNC2(v(mV)){
	UNITSOFF
      tmUNC2=(p1tmunc2/(exp(-(v-p2tmunc2+shiftmunc2)/p3tmunc2)+exp((v-p2tmunc2+shiftmunc2)/p4tmunc2))+p5tmunc2)*constmunc2
     UNITSON
	  }

FUNCTION thunc2(v(mV)){
     UNITSOFF
      thunc2=(p1thunc2/(1+exp((v-p2thunc2+shifthunc2)/(p3thunc2*fp5)))+p4thunc2/(1+exp(-(v-p5thunc2+shifthunc2)/(p6thunc2*fp5))))*consthunc2
     UNITSON
	  }
