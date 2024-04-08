TITLE caintra1
: Intracellular calcium for sk channels
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)
: "Biophysical modeling of the whole-cell dynamics of C. elegans motor and interneurons families"
: M. Nicoletti et al. PloS ONE, 19(3): e0298105.
: https://doi.org/10.1371/journal.pone.0298105


NEURON {
SUFFIX caintra1
USEION ca READ ica,eca 
RANGE vol, surf
GLOBAL ca_eq,Fc,calcium
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	PI = (pi) (1)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um)= (micrometer)

}




PARAMETER {
	 vol=7.42e-12 (cm3)
	 surf=65.89e-8 (cm2)
	 fca = 0.001
	 tca = 50 (ms) 
	 ca_eq= 0.05e-6 (M)
	 eca (mV)
	 Fc=96485 (coul)

}


ASSIGNED { ica (mA/cm2) 
        calcium (M)
 
}


STATE { caintra }

INITIAL{

caintra=ca_eq

}

BREAKPOINT {
SOLVE state METHOD cnexp 
    calcium=caintra
}

DERIVATIVE state {
LOCAL rs
if (ica<=0){

rs=fca*(-((1/(2*vol*Fc))*(ica*surf*1e-3)))-((caintra-ca_eq)/tca)

} else {
rs=-((caintra-ca_eq)/tca)
}

caintra'=rs
}

