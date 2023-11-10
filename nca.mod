TITLE NCA
: A passive leak current
: From Nicoletti et al. PloS One 2019 (https://doi.org/10.1371/journal.pone.0218738)

UNITS {
    (mA) = (milliamp)
	(S) = (siemens)
	(mV) = (millivolt)
}

NEURON {
	SUFFIX nca
	NONSPECIFIC_CURRENT i
	RANGE i, e, gbar
}



PARAMETER {
	gbar = 0.055  (S/cm2)  
	e = 30  (mV)
}

ASSIGNED {
	i  (mA/cm2)
	v  (mV)
}

BREAKPOINT { i = gbar*(v - e) }