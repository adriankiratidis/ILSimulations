!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program run_C4MIM_BF4
  use ILsimulationssrclib
  implicit none

  character(len=256) :: file_stub

  print *, "Please Enter the data file prefix"
  read(*,*) file_stub
  
  print *, "Reading and initialising model parameters"
  print *, "This includes discretisation params and simulation params"
  call InitialiseModelParameters(trim(file_stub))

  
  

end program run_C4MIM_BF4



  
