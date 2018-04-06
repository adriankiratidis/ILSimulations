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
  ! over range of hs_diameter. (i.e. >= x and < x + hs_diameter)
  public :: n_discretised_points_z

  public :: epsilonr
  public :: epsilon_LJ
  public :: surface_charge_density
  public :: hs_diameter
  public :: plate_separations ! in multiples of hs_diameter
  public :: bulk_density
  public :: bead_charge
  public :: string_length
  public :: iterative_tolerance
  public :: max_iteration_limit

  real(dp) :: chi_parameter
  integer :: n_discretised_points_z

  real(dp) :: epsilonr
  real(dp) :: epsilon_LJ
  real(dp) :: surface_charge_density
  real(dp) :: hs_diameter

  !array of plater separations in multiples of hs_diameter
  integer, dimension(:), allocatable  :: plate_separations
  
  real(dp) :: bulk_density
  real(dp) :: bead_charge
  real(dp) :: string_length
  real(dp) :: iterative_tolerance
  integer  :: max_iteration_limit


contains

  !Subroutine that reads in all required params
  subroutine InitialiseModelParameters(file_stub)
    character(len=*), intent(in) :: file_stub

    integer :: iseparation
    integer :: n_plate_separations
    integer :: file_unit
    file_unit = 171

    open(file_unit, file=trim(file_stub)//".params", action='read')

    read(file_unit, *) chi_parameter
    read(file_unit, *) epsilonr
    read(file_unit, *) surface_charge_density
    read(file_unit, *) hs_diameter
    read(file_unit, *) bulk_density
    read(file_unit, *) bead_charge
    read(file_unit, *) string_length
    read(file_unit, *) n_discretised_points_z
    read(file_unit, *) max_iteration_limit
    read(file_unit, *) iterative_tolerance
    read(file_unit, *) n_plate_separations

    allocate(plate_separations(n_plate_separations))
    do iseparation = 1, n_plate_separations
       read(file_unit, *) plate_separations(iseparation)
    end do

    close(file_unit)

    print *,  "Succesfully set the following values"
    print *,  "chi_parameter = ", chi_parameter
    print *,  "epsilonr = ", epsilonr
    print *,  "surface_charge_density = ", surface_charge_density
    print *,  "hs_diameter = ", hs_diameter
    print *,  "bulk_density = ", bulk_density
    print *,  "bead_charge = ", bead_charge
    print *,  "string_length = ", string_length
    print *,  "n_discretised_points_z = ", n_discretised_points_z
    print *,  "max_iteration_limit = ", max_iteration_limit
    print *,  "iterative tolerance = ", iterative_tolerance
    print *,  "n_plate_separations = ", n_plate_separations
    print *,  "plate separations are: ", plate_separations

  end subroutine InitialiseModelParameters

  subroutine DeAllocateModelParams()

    if(allocated(plate_separations)) deallocate(plate_separations)

  end subroutine DeAllocateModelParams
end module parameters
