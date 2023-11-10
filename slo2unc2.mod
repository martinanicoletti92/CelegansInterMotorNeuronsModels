TITLE slo2unc2
: slo2 channels coupled with unc2 calcium channels (1:1 stoichiometry)
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738) 

UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
	 (molar)=(1/liter)
	  FARADAY = (faraday) (coulombs)
	 (uM) =	(micromolar)
	 (pS) = (picosiemens)
}

NEURON {
	SUFFIX slo2unc2
	USEION k READ ek WRITE ik
	USEION ca READ eca
	RANGE  gbar
	EXTERNAL  munc2_unc2, hunc2_unc2
 }


PARAMETER{
   v (mV)
   cai (uM)
   fondo=0.05 (uM)
   ek (mV)
   eca (mV)

   munc2_unc2
   hunc2_unc2
  
   celsius		(degC)
   gbar=0.1 (S/cm2)
   wom1=0.896395  (/ms)  
   wyx1=0.019405  (/mV)
   kyx1=3294.553404 (uM)
   nyx1=0.000010   (1)
   wop1=0.026719 (/ms)
   wxy1=-0.024123 (/mV)
   kxy1=93.449423 (/ms)
   nxy1=1.835067 (1)
   
   r=13e-9 (m)
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
		 minf
		 tslo2
		 alpha
         beta		 
}	


STATE {
	m 
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*m*hunc2_unc2*(v-ek)
}


INITIAL {
	rates(calcio(v), v)
	m=minf
}

DERIVATIVE states { 
       rates(calcio(v), v)   
        m' = (minf - m)/tslo2
	
	}


PROCEDURE rates(calcio(v),v (mV)){ 
    alpha=minfunc2(v)/tmunc2(v)
	beta=(1/tmunc2(v))-alpha
 	minf=(munc2_unc2*kop2(calcio(v),v)*(alpha+beta+kcm2(v)))/((kop2(calcio(v),v)+kom2(calcio(v),v))*(kcm2(v)+alpha)+(beta*kcm2(v)))
	tslo2=((alpha+beta+kcm2(v))/((kop2(calcio(v),v)+kom2(calcio(v),v))*(kcm2(v)+alpha)+(beta*kcm2(v))))
 }

		
FUNCTION kcm2(v (mV)){
		kcm2=wom1*exp(-wyx1*v)*(1/(1+pow(fondo/kyx1,nyx1)))
}
	
FUNCTION kom2(calcio(v),v (mV)){
		kom2=wom1*exp(-wyx1*v)*(1/(1+pow(calcio(v)/kyx1,nyx1)))
}	

FUNCTION kop2(calcio(v),v (mV)){
		kop2=wop1*exp(-wxy1*v)*(1/(1+pow(kxy1/calcio(v),nxy1)))
}



FUNCTION calcio(v (mV)){
	calcio=(((fabs(gsc*(v-eca)*1e-3)/(8*pi*r*d*FARADAY))*exp(-r/sqrt(d/(kb*b))))*1e6*1e-3)+fondo
}

FUNCTION minfunc2(v(mV)){

     minfunc2=1/(1+exp(-(v-va_unc2+stm2)/(ka_unc2*func2)))
}

FUNCTION hinfunc2(v(mV)){

     hinfunc2= 1/(1+exp((v-vi_unc2+sth2)/(ki_unc2*f2unc2)))
  
}

FUNCTION tmunc2(v(mV)){

      tmunc2=(p1tmunc2/(exp(-(v-p2tmunc2+shiftmunc2)/(p3tmunc2*fp3))+exp((v-p2tmunc2+shiftmunc2)/(p4tmunc2*fp4)))+p5tmunc2)*constmunc2
     
	  }

FUNCTION thunc2(v(mV)){
     
      thunc2=(p1thunc2/(1+exp((v-p2thunc2+shifthunc2)/(p3thunc2*fp5)))+p4thunc2/(1+exp(-(v-p5thunc2+shifthunc2)/(p6thunc2*fp5))))*consthunc2
    
	  }

