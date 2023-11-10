TITLE slo2egl19
: slo2 channels coupled with egl19 calcium channels (1:1 stoichiometry)
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738) 

UNITS {
(mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
     (pS) = (picosiemens)
	 (molar)=(1/liter)
	 (uM) =	(micromolar)
	 (um) = (micron)
	 FARADAY = (faraday) (coulombs)
}

NEURON {
	SUFFIX slo2egl19
	USEION k READ ek WRITE ik
	USEION ca READ eca
	RANGE  gbar, g, curr
	GLOBAL minf,tslo2, alpha1,beta1
	EXTERNAL  megl19_egl19, hegl19_egl19
 }


PARAMETER{
   v (mV)
   cai (uM)
   bkg=0.05 (uM)
   ek (mV)
   eca (mV)
   minf_egl19
   hinf_egl19
   mtau_egl192
   htau_egl192
   megl19_egl19
   hegl19_egl19
  
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
   
   r=13e-9 (nm)
   d=250e-12 (um2/s)
   kb=500e6 (1/(M-s)
   b=30e-6 (M)
   gsc=40e-12 (pS)
   pi=3.14
   
   va_egl19=5.6 (mV)
	ka_egl19=7.50 (mV)
	stm19=10 (mV)
	sth19=10 (mV)
   
    p1hegl19=1.4314 
    p2hegl19=24.8573 (mV)
    p3hegl19=11.9541 (mV)
    p4hegl19=0.1427
    p5hegl19=5.9589
    p6hegl19=-10.5428 (mV)
    p7hegl19=8.0552 (mV)
    p8hegl19=0.6038
  
    pdg1=2.3359 (ms)
	pdg2=2.9324 (ms)
	pdg3=5.2357 (mV)
	pdg4=6.0 (mV)
	pdg5=1.8739 (ms)
	pdg6=1.3930 (mV)
	pdg7=30.0 (mV)
	stau19=10 (mV)  

  
	
	
	pds1=0.4 
	pds2=0.55 
	pds3=81.1179  (ms)
	pds4=-22.9723 (mV)
	pds5=5 (mV)
	pds6=43.0937 (ms)
	pds7=0.9
	pds8=40.4885 (ms)
	pds9=28.7251 (mV)
	pds10=3.7125 (mV)
	pds11=0
	shiftdps=10
	ctm19=1
   
    }

ASSIGNED{
		 ik (mA/cm2)
		 g (S/cm2)
		 curr (mA/cm2)
		 minf
		 tslo2
		 alpha1
         beta1		 
}	


STATE {
	m 
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g=gbar*m*hegl19_egl19
	curr=gbar*m*hegl19_egl19*(v-ek)
	ik = gbar*m*hegl19_egl19*(v-ek)
}


INITIAL {
	rates(calcium(v), v)
	m=minf
}

DERIVATIVE states { 
       rates(calcium(v), v)   
        m' = (minf - m)/tslo2
	
	}


PROCEDURE rates(calcium(v),v (mV)){ 
    alpha1=actegl19(v)/tactegl19(v)
	beta1=(1/tactegl19(v))-alpha1
 	minf=(megl19_egl19*kop2(calcium(v),v)*(alpha1+beta1+kcm2(v)))/((kop2(calcium(v),v)+kom2(calcium(v),v))*(kcm2(v)+alpha1)+(beta1*kcm2(v)))
	tslo2=((alpha1+beta1+kcm2(v))/((kop2(calcium(v),v)+kom2(calcium(v),v))*(kcm2(v)+alpha1)+(beta1*kcm2(v))))
 }

		
FUNCTION kcm2(v (mV)){
		kcm2=wom1*exp(-wyx1*v)*(1/(1+((bkg/kyx1)^nyx1)))
}
	
FUNCTION kom2(calcium(v),v (mV)){
		kom2=wom1*exp(-wyx1*v)*(1/(1+pow(calcium(v)/kyx1,nyx1)))
}	

FUNCTION kop2(calcium(v),v (mV)){
		kop2=wop1*exp(-wxy1*v)*(1/(1+pow(kxy1/calcium(v),nxy1)))
}

FUNCTION calcium(v (mV)){
	calcium=(((fabs(gsc*(v-eca)*1e-3)/(8*pi*r*d*FARADAY))*exp(-r/sqrt(d/(kb*b))))*1e6*1e-3)+bkg
	
}


FUNCTION actegl19(v (mV)) { 
		UNITSOFF
        actegl19=1/(1+exp(-(v-va_egl19+stm19)/ka_egl19))
		UNITSON
		}
FUNCTION inactegl19(m (mV)){
		UNITSOFF
        inactegl19= ((p1hegl19/(1+exp(-(v-p2hegl19+sth19)/p3hegl19))+p4hegl19)*(p5hegl19/(1+exp((v-p6hegl19+sth19)/p7hegl19))+p8hegl19))
        UNITSON
		}
FUNCTION tactegl19(m (mV)){
		UNITSOFF
		tactegl19= (pdg1+(pdg2*exp(-(v-pdg3+stau19)^2/(pdg4)^2))+(pdg5*exp(-(v-pdg6+stau19)^2/(pdg7)^2)))*ctm19
		UNITSON
		}
		
FUNCTION tinactegl19(v(mV)){
		UNITSOFF
		tinactegl19= pds1*(((pds2*pds3)/(1+exp((v-pds4+shiftdps)/pds5)))+pds6+((pds7*pds8)/(1+exp((v-pds9+shiftdps)/pds10)))+pds11)
		UNITSON
		}


