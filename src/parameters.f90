!Module storing routines that read and initialise all required parameters.
module parameters
  use kinds
  implicit none
  public

  public :: InitialiseModelParameters

  ! 'd' = chi_parameter * hs_diameter provides the only adjusted
  ! parameter in the model. It is fitted to give the correct bulk
  ! ion density for the particlar RTIL at a given pressure and temp.
  ! Only present in the unlike electric correlation term.
  public :: chi_parameter

  ! Number of discretised points in the z direction
  ! over range of hs_diameter. (i.e. between x and x + hs_diameter)
  public :: n_discretised_points_z

  public :: epsilonr
  public :: surface_charge_density
  public :: hs_diameter
  public :: plate_separation
  public :: bulk_density
  public :: string_length

  real(dp) :: chi_parameter
  integer :: n_discretised_points_z
  
  real(dp) :: epsilonr
  real(dp) :: surface_charge_density
  real(dp) :: hs_diameter
  integer  :: plate_separation !In multiples of hs_diameter
  real(dp) :: bulk_density
  real(dp) :: string_length
  integer  :: MAX_ITERATION_LIMIT
  
contains

  !Subroutine that reads in all required params
  subroutine InitialiseModelParameters(file_stub)
    character(len=*), intent(in) :: file_stub

    integer :: file_unit
    file_unit = 171

    open(file_unit, file=trim(file_stub)//".params", action='read')
    read(file_unit, *) chi_parameter
    read(file_unit, *) epsilonr
    read(file_unit, *) surface_charge_density
    read(file_unit, *) hs_diameter
    read(file_unit, *) plate_separation
    read(file_unit, *) bulk_density
    read(file_unit, *) string_length
    read(file_unit, *) n_discretised_points_z
    read(file_unit, *) MAX_ITERATION_LIMIT
    close(file_unit)

    print *, "Succesfully set the following values"
    print *,  "chi_parameter = ", chi_parameter
    print *,  "epsilonr = ", epsilonr
    print *,  "surface_charge_density = ", surface_charge_density
    print *,  "hs_diameter = ", hs_diameter
    print *,  "plate_separation = ", plate_separation
    print *,  "bulk_density = ", bulk_density
    print *,  "string_length = ", string_length
    print *,  "n_discretised_points_z = ", n_discretised_points_z
    print *, "MAX_ITERATION_LIMIT = ", MAX_ITERATION_LIMIT

  end subroutine InitialiseModelParameters

end module parameters
