TITLE slo1egl19
: slo1 channels coupled with egl19 calcium channels (1:1 stoichiometry)
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)

UNITS {
	 (mA) = (milliamp)
	 (mV) = (millivolt)
	 (nS) = (nanosiemens)
     (pS) = (picosiemens)
	 (molar)=(1/liter)
	 (uM) =	(micromolar)
	  FARADAY = (faraday) (coulombs)
	 (M)= (molar)
}

NEURON {
	SUFFIX slo1egl19
	USEION k READ ek WRITE ik
	USEION ca READ eca
	GLOBAL alpha1,beta1
	RANGE  gbar,g,mminf,tslo1,alph,bet,p,mca,hca,taca,tica,mminfca,hinfca,kcmCALC,kopCALC,komCALC,caCALC,m1,m2,curr
	EXTERNAL  megl19_egl19, hegl19_egl19
 }


PARAMETER{
   v (mV)
   cai (uM)
   fondo=0.05 (uM)
   ek (mV)
   eca (mV)
   megl19_egl19
   hegl19_egl19
  
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
   
   r=13e-9 (m)
   d=250e-12 (um2/s)
   kb=500e6 (/M-s)
   b=30e-6 (M)
   gsc=40e-12 (S)
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
		 curr(mA/cm2)
		 mminf
		 tslo1
		 alph
		 alpha1
		 beta1
         bet
		 mca
		 hca
		 taca
		 tica
		 mminfca
		 hinfca
		 kcmCALC
		 kopCALC
		 komCALC
		 caCALC
		 m1
		 m2
		 p
}	


STATE {
	m 
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g=gbar*m*hegl19_egl19
	mminf=mminf
	tslo1=tslo1
	mca=megl19_egl19
	hca=hegl19_egl19
	taca=tactegl19(v)
	tica=tinactegl19(v)
	mminfca=actegl19(v)
	hinfca=inactegl19(v)
    kcmCALC=kcm(v)
    komCALC=kom(calcium(v),v)
    kopCALC=kop(calcium(v),v)
	caCALC=calcium(v)
	m1=megl19_egl19
	m2=hegl19_egl19
	alph=alpha1
	bet=beta1

	curr=gbar*m*hegl19_egl19*(v-ek)
	ik = gbar*m*hegl19_egl19*(v-ek)
}


INITIAL {
	rates(calcium(v), v)
	m=mminf
}

DERIVATIVE states { 
       rates(calcium(v), v)   
        m' = (mminf - m)/tslo1
	
	}


PROCEDURE rates(calcium(v),v (mV)){ 
    alpha1=actegl19(v)/tactegl19(v)
	beta1=(1/tactegl19(v))-alpha1
 	mminf=(megl19_egl19*kop(calcium(v),v)*(alpha1+beta1+kcm(v)))/((kop(calcium(v),v)+kom(calcium(v),v))*(kcm(v)+alpha1)+(beta1*kcm(v)))
	tslo1=((alpha1+beta1+kcm(v))/((kop(calcium(v),v)+kom(calcium(v),v))*(kcm(v)+alpha1)+(beta1*kcm(v))))
 }

		
FUNCTION kcm(v (mV)){
		kcm=wom*exp(-wyx*v)*(1/(1+((fondo/kyx)^nyx)))
}
	
FUNCTION kom(calcium(v),v (mV)){
		kom=wom*exp(-wyx*v)*(1/(1+pow(calcium(v)/kyx,nyx)))
}	

FUNCTION kop(calcium(v),v (mV)){
		kop=wop*exp(-wxy*v)*(1/(1+pow(kxy/calcium(v),nxy)))
}

FUNCTION calcium(v (mV)){
	calcium=(((fabs(gsc*(v-eca)*1e-3)/(8*pi*r*d*FARADAY))*exp(-r/sqrt(d/(kb*b))))*1e6*1e-3)+fondo
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
