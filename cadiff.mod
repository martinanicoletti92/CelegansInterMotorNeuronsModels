: Ca diffusion in a Purkinje cell
: Created 8/15/02 - nwg
: https://senselab.med.yale.edu/ModelDB/ShowModel?model=48332&file=/purkinje/cadiff.mod#tabs-1


NEURON {
       SUFFIX cadiff
       USEION ca READ ica, cai WRITE cai
       RANGE ca
       GLOBAL depth, beta
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
      (um) = (micron)
}

CONSTANT {
      F = 9.6485e4 (coul)
}

PARAMETER {
          cai      (mM)
          dt       (ms)

          depth  = .1  (um)
          beta = 1 (/ms)
}

ASSIGNED {
         ica       (mA/cm2)
}

STATE {
      ca           (mM)
}

INITIAL {
        ca = .0001
}

BREAKPOINT {

        ca = ca + (10000) * dt * ( ( -1/(2*F)*ica / (depth)) - (.0001) * beta * ca )

        if ( ca < 1e-4 ) {: minimum 100 nM Ca
           ca = 1e-4
        }

        cai = ca

}