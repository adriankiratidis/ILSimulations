!Module storing routines that read and initialise all required parameters.
module parameters
  use kinds
  use universalconstants
  implicit none
  public

  !Public Subroutines
  public :: InitialiseModelParameters
  public :: SetBulkDensityOfEndMonomers

  !Public Variables
  public :: ionic_liquid_name 

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
  public :: surface_charge_density_left_wall
  public :: surface_charge_density_right_wall
  public :: hs_diameter
  public :: a_term_index !Index on a_ex in hard sphere term
  public :: plate_separations ! in multiples of hs_diameter
  public :: bulk_density
  public :: temperature
  public :: positive_bead_charge
  public :: negative_bead_charge
  public :: string_length
  public :: iterative_tolerance
  public :: max_iteration_limit
  public :: beta
  public :: alpha_mixing_for_update
  
  public :: bulk_density_positive_beads
  public :: bulk_density_neutral_beads
  public :: bulk_density_negative_beads
  
  public :: r_plus
  public :: r_neutral
  public :: r_minus
  
  character(len=256) :: ionic_liquid_name 
  real(dp) :: chi_parameter
  integer  :: n_discretised_points_z
  real(dp) :: epsilonr
  real(dp) :: epsilon_LJ
  real(dp) :: surface_charge_density_left_wall
  real(dp) :: surface_charge_density_right_wall
  real(dp) :: hs_diameter
  integer  :: a_term_index

  !array of plater separations in multiples of hs_diameter
  real(dp), dimension(:), allocatable  :: plate_separations

  real(dp) :: bulk_density
  real(dp) :: temperature
  real(dp) :: positive_bead_charge
  real(dp) :: negative_bead_charge
  real(dp) :: string_length
  real(dp) :: iterative_tolerance
  integer  :: max_iteration_limit
  real(dp) :: beta
  real(dp) :: alpha_mixing_for_update
  
  integer :: r_plus
  integer :: r_neutral
  integer :: r_minus
  
  real(dp) :: bulk_density_positive_beads
  real(dp) :: bulk_density_neutral_beads
  real(dp) :: bulk_density_negative_beads

contains

  !Subroutine that reads in all required params
  subroutine InitialiseModelParameters(file_stub)
    character(len=*), intent(in) :: file_stub

    integer :: iseparation
    integer :: n_plate_separations
    integer :: file_unit
    file_unit = 171

    open(file_unit, file=trim(file_stub)//".params", action='read')

    read(file_unit, *) ionic_liquid_name
    read(file_unit, *) chi_parameter
    read(file_unit, *) epsilonr
    read(file_unit, *) epsilon_LJ !
    read(file_unit, *) surface_charge_density_left_wall
    read(file_unit, *) surface_charge_density_right_wall
    read(file_unit, *) hs_diameter
    read(file_unit, *) a_term_index
    read(file_unit, *) bulk_density !
    read(file_unit, *) temperature
    read(file_unit, *) alpha_mixing_for_update
    read(file_unit, *) positive_bead_charge
    read(file_unit, *) negative_bead_charge
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

    ! Set derived parameters
    beta = 1.0_dp / (k_B * temperature)

    ! Apply unit transformations
    epsilon_LJ = epsilon_LJ !* k_B
    bulk_density = bulk_density / (hs_diameter**3.0_dp)

    print *,  "Succesfully set the following values"
    print *,  "ionic_liquid_name = ", ionic_liquid_name
    print *,  "chi_parameter = ", chi_parameter
    print *,  "epsilonr = ", epsilonr
    print *,  "epsilon_LJ = ", epsilon_LJ
    print *,  "surface_charge_density_left_wall = ", surface_charge_density_left_wall
    print *,  "surface_charge_density_right_wall = ", surface_charge_density_right_wall
    print *,  "hs_diameter = ", hs_diameter
    print *,  "a_term_index = ", a_term_index
    print *,  "bulk_density = ", bulk_density
    print *,  "temperature = ", temperature
    print *,  "alpha_mixing_for_update = ", alpha_mixing_for_update
    print *,  "positive_bead_charge = ", positive_bead_charge
    print *,  "negative_bead_charge = ", negative_bead_charge
    print *,  "string_length = ", string_length
    print *,  "n_discretised_points_z = ", n_discretised_points_z
    print *,  "max_iteration_limit = ", max_iteration_limit
    print *,  "iterative tolerance = ", iterative_tolerance
    print *,  "n_plate_separations = ", n_plate_separations
    print *,  "plate separations are: ", plate_separations
    print *,  "thermodynamic beta = ", beta

    print *, "Setting Bead Densities from the Bulk Ion Density."
    call SetBeadDensityFromBulkIonDensity(ionic_liquid_name)

    print *, "Checking plate separations are valid given discretisation scheme."
    call CheckValidityOfPlateSeparations()

  end subroutine InitialiseModelParameters

  subroutine DeAllocateModelParams()

    if(allocated(plate_separations)) deallocate(plate_separations)

  end subroutine DeAllocateModelParams

  subroutine SetBeadDensityFromBulkIonDensity(ionic_liquid_name)
    character(len=*), intent(in) :: ionic_liquid_name

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       call SetSingleNeutralSphereBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then
       call SetSinglePositiveSphereBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then
       call SetSingleNegativeSphereBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       call SetNeutralDimerDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then
       call SetC4MIN_BF4BeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then
       call SetPositiveMinusBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then
       call SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then
       call SetDimerDoubleDimerBeadDensityFromBulkIonDensity()


    else
       print *, "parameters.f90: SetBeadDensityFromBulkIonDensity:"
       print *, "Unsupported 'ionic_liquid_name' value of ", trim(ionic_liquid_name)
       print *, "...aborting..."
       call abort()
    end if

  end subroutine SetBeadDensityFromBulkIonDensity

  subroutine SetSingleNeutralSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = 0.0_dp

    r_plus = 0
    r_neutral = 0
    r_minus = 0

  end subroutine SetSingleNeutralSphereBeadDensityFromBulkIonDensity

  subroutine SetSinglePositiveSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = 0.0_dp
  end subroutine SetSinglePositiveSphereBeadDensityFromBulkIonDensity

  subroutine SetSingleNegativeSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = bulk_density
  end subroutine SetSingleNegativeSphereBeadDensityFromBulkIonDensity

  subroutine SetNeutralDimerDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = 2.0_dp * bulk_density
    bulk_density_negative_beads = 0.0_dp
  end subroutine SetNeutralDimerDensityFromBulkIonDensity

  subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 5.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density

    r_plus = 3
    r_neutral = 3
    r_minus = 1

  end subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetPositiveMinusBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = bulk_density
  end subroutine SetPositiveMinusBeadDensityFromBulkIonDensity

  subroutine SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = bulk_density

  end subroutine SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity

  subroutine SetDimerDoubleDimerBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = 2.0_dp * bulk_density
    bulk_density_neutral_beads = 2.0_dp * bulk_density
    bulk_density_negative_beads = 2.0_dp * bulk_density

    r_plus = 1
    r_neutral = 1
    r_minus = 0
  end subroutine SetDimerDoubleDimerBeadDensityFromBulkIonDensity

  subroutine SetBulkDensityOfEndMonomers(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, n_minus_total, n_minus_end)
    real(dp), dimension(:), intent(in) :: n_plus_total
    real(dp), dimension(:), intent(out) :: n_plus_end
    real(dp), dimension(:), intent(in) :: n_neutral_total
    real(dp), dimension(:), intent(out) :: n_neutral_end
    real(dp), dimension(:), intent(in) :: n_minus_total
    real(dp), dimension(:), intent(out) :: n_minus_end

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       n_plus_end = n_plus_total
       n_neutral_end = n_neutral_total
       n_minus_end = n_minus_total

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then
       n_plus_end = n_plus_total
       n_neutral_end = n_neutral_total
       n_minus_end = n_minus_total

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then
       n_plus_end = n_plus_total
       n_neutral_end = n_neutral_total
       n_minus_end = n_minus_total

    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       n_plus_end = n_plus_total
       n_neutral_end = n_neutral_total
       n_minus_end = n_minus_total

    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then
       n_plus_end = 0.4_dp * n_plus_total
       n_neutral_end = 0.4_dp * n_neutral_total
       n_minus_end = 0.8_dp * n_minus_total

    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then
       n_plus_end = n_plus_total
       n_neutral_end = n_neutral_total
       n_minus_end = n_minus_total

    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then
       n_plus_end = n_plus_total
       n_neutral_end = n_neutral_total
       n_minus_end = n_minus_total

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then
       n_plus_end = 0.5_dp * n_plus_total
       n_neutral_end = 0.5_dp * n_neutral_total
       n_minus_end = n_minus_total

    else
       print *, "parameters.f90: SetBeadDensityFromBulkIonDensity:"
       print *, "Unsupported 'ionic_liquid_name' value of ", trim(ionic_liquid_name)
       print *, "...aborting..."
       call abort()
    end if

  end subroutine SetBulkDensityOfEndMonomers



  subroutine CheckValidityOfPlateSeparations()

    integer :: ith_separation

    real(dp) :: distance_between_z_values
    real(dp) :: distance_beyond_hs_diameter_multiple
    real(dp) :: n_extra_discretisation_points

    distance_between_z_values = hs_diameter / n_discretised_points_z

    !Check that all the separations are a multiple of the distance between allowed points
    do ith_separation = 1, size(plate_separations)

       distance_beyond_hs_diameter_multiple = hs_diameter * &
            (plate_separations(ith_separation) - floor(plate_separations(ith_separation)))

       n_extra_discretisation_points = distance_beyond_hs_diameter_multiple / distance_between_z_values

       if(abs(n_extra_discretisation_points - nint(n_extra_discretisation_points)) > 1.0E-6) then
          print *, "parameter.f90: CheckValidityOfPlateSeparations: "
          print *, "abs(n_extra_discretisation_points - nint(n_extra_discretisation_points)) > 1.0E-6"
          print *, "abs(n_extra_discretisation_points - nint(n_extra_discretisation_points)) = ", &
               abs(n_extra_discretisation_points - nint(n_extra_discretisation_points))
          print *, "distance_beyond_hs_diameter_multiple = ", distance_beyond_hs_diameter_multiple
          print *, "n_extra_discretisation_points = ", n_extra_discretisation_points
          print *, "Total distance between plates is not an integer multiple of the distance between"
          print *, "discretisations.  Input error.  Change plate separations...aborting..."
          call abort()
       end if

    end do

  end subroutine CheckValidityOfPlateSeparations

end module parameters
