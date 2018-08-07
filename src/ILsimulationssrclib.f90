!Module that acts as a library that stores all the source modules.
!To be used by programs requiring the entire src.
module ILsimulationssrclib
  use kinds
  use universalconstants
  use parameters
  use integratephispherical
  use integratezcylindrical
  use iteration
  use lambdas
  use helpers
  use io
  use normalisation
  use surfaceforces
  use constructoligomers
  use contacttheorem
  use discretederivatives
  use excessenergyfunctionalparameters
  use chargeincrements
  implicit none
end module ILsimulationssrclib
