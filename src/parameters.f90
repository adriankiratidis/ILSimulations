!Module storing routines that read and initialise all required parameters.
module parameters
  use kinds
  use universalconstants
  implicit none
  public

  public :: InitialiseModelParameters

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
  public :: epsilon_LJ_particle_particle
  public :: epsilon_LJ_particle_wall
  public :: epsilon_eighth_power_const
  public :: type_of_interaction_r4r8
  public :: mica_density
  public :: surface_charge_density_left_wall
  public :: surface_charge_density_right_wall
  public :: hs_diameter
  public :: a_term_index !Index on a_ex in hard sphere term
  public :: plate_separations ! in multiples of hs_diameter
  public :: bulk_density
  public :: temperature
  public :: positive_bead_charge
  public :: negative_bead_charge
  public :: positive_oligomer_charge
  public :: negative_oligomer_charge
  public :: string_length
  public :: iterative_tolerance
  public :: max_iteration_limit
  public :: beta
  public :: alpha_mixing_for_update

  public :: bulk_density_positive_beads
  public :: bulk_density_neutral_beads
  public :: bulk_density_negative_beads

  public :: slope_for_initial_guess
  public :: n_charge_iterations
  public :: Donnan_potential

  public :: r_cation !r parameter used in hard sphere paper in calculating Ymix
  public :: r_anion
  public :: num_end_monomers_cation
  public :: num_end_monomers_anion

  public :: n_plus_cation_end_bulk
  public :: n_neutral_cation_end_bulk
  public :: n_minus_cation_end_bulk

  public :: n_plus_cation_nonend_bulk
  public :: n_neutral_cation_nonend_bulk
  public :: n_minus_cation_nonend_bulk

  public :: n_plus_anion_end_bulk
  public :: n_neutral_anion_end_bulk
  public :: n_minus_anion_end_bulk

  character(len=256) :: ionic_liquid_name 
  real(dp) :: chi_parameter
  integer  :: n_discretised_points_z
  real(dp) :: epsilonr
  real(dp) :: epsilon_LJ_particle_particle
  real(dp) :: epsilon_LJ_particle_wall
  real(dp) :: epsilon_eighth_power_const
  character(len=16) :: type_of_interaction_r4r8
  real(dp) :: mica_density
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
  real(dp) :: positive_oligomer_charge
  real(dp) :: negative_oligomer_charge
  real(dp) :: string_length
  real(dp) :: iterative_tolerance
  integer  :: max_iteration_limit
  real(dp) :: beta
  real(dp) :: alpha_mixing_for_update

  real(dp) :: bulk_density_positive_beads
  real(dp) :: bulk_density_neutral_beads
  real(dp) :: bulk_density_negative_beads

  real(dp) :: slope_for_initial_guess
  integer :: n_charge_iterations

  real(dp) :: Donnan_potential

  integer :: r_cation
  integer :: r_anion
  integer :: num_end_monomers_cation
  integer :: num_end_monomers_anion

  real(dp) :: n_plus_cation_end_bulk
  real(dp) :: n_neutral_cation_end_bulk
  real(dp) :: n_minus_cation_end_bulk

  real(dp) :: n_plus_cation_nonend_bulk
  real(dp) :: n_neutral_cation_nonend_bulk
  real(dp) :: n_minus_cation_nonend_bulk

  real(dp) :: n_plus_anion_end_bulk
  real(dp) :: n_neutral_anion_end_bulk
  real(dp) :: n_minus_anion_end_bulk

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
    read(file_unit, *) epsilon_LJ_particle_particle !
    read(file_unit, *) epsilon_LJ_particle_wall !
    read(file_unit, *) epsilon_eighth_power_const
    read(file_unit, *) type_of_interaction_r4r8
    read(file_unit, *) mica_density
    read(file_unit, *) surface_charge_density_left_wall
    read(file_unit, *) surface_charge_density_right_wall
    read(file_unit, *) hs_diameter
    read(file_unit, *) a_term_index
    read(file_unit, *) bulk_density !
    read(file_unit, *) temperature
    read(file_unit, *) alpha_mixing_for_update
    read(file_unit, *) slope_for_initial_guess
    read(file_unit, *) n_charge_iterations
    read(file_unit, *) positive_bead_charge
    read(file_unit, *) negative_bead_charge
    read(file_unit, *) Donnan_potential !Initial Guess for the Donnan potential
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
    epsilon_LJ_particle_particle = (epsilon_LJ_particle_particle * k_B) !/ ((hs_diameter**3)*mica_density)
    epsilon_LJ_particle_wall = epsilon_LJ_particle_wall * k_B * ((hs_diameter**3)*mica_density)

    bulk_density = bulk_density / (hs_diameter**3.0_dp)

    positive_bead_charge = positive_bead_charge * electric_charge
    negative_bead_charge = negative_bead_charge * electric_charge

    surface_charge_density_left_wall = surface_charge_density_left_wall * electric_charge
    surface_charge_density_right_wall = surface_charge_density_right_wall * electric_charge


    print *,  "Succesfully set the following values"
    print *,  "ionic_liquid_name = ", ionic_liquid_name
    print *,  "chi_parameter = ", chi_parameter
    print *,  "epsilonr = ", epsilonr
    print *,  "epsilon_LJ particle_particle interaction = ", epsilon_LJ_particle_particle
    print *,  "epsilon_LJ particle_wall interaction = ", epsilon_LJ_particle_wall
    print *,  "epsilon_eighth_power_const = ", epsilon_eighth_power_const
    print *,  "mica density = ", mica_density
    print *,  "surface_charge_density_left_wall = ", surface_charge_density_left_wall
    print *,  "surface_charge_density_right_wall = ", surface_charge_density_right_wall
    print *,  "hs_diameter = ", hs_diameter
    print *,  "a_term_index = ", a_term_index
    print *,  "bulk_density = ", bulk_density
    print *,  "temperature = ", temperature
    print *,  "alpha_mixing_for_update = ", alpha_mixing_for_update
    print *,  "slope_for_initial_guess = ", slope_for_initial_guess
    print *,  "n_charge_iterations = ", n_charge_iterations
    print *,  "positive_bead_charge = ", positive_bead_charge
    print *,  "negative_bead_charge = ", negative_bead_charge
    print *,  "Initial Guess for the Donnan potential = ", Donnan_potential
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

    else if(trim(ionic_liquid_name) == "C2MIM_BF4-") then
       call SetC2MIN_BF4BeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C6MIM_BF4-") then
       call SetC6MIN_BF4BeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C8MIM_BF4-") then
       call SetC8MIN_BF4BeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C10MIM_BF4-") then
       call SetC10MIN_BF4BeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C4MIM+_TFSI-_model1") then
       call SetC4MIN_TFSIBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C2MIM+_TFSI-_model1") then
       call SetC2MIN_TFSIBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C6MIM+_TFSI-_model1") then
       call SetC6MIN_TFSIBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C8MIM+_TFSI-_model1") then
       call SetC8MIN_TFSIBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C10MIM+_TFSI-_model1") then
       call SetC10MIN_TFSIBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C4MIM+_TFSI-_model2") then
       call SetC4MIN_TFSIBeadDensityFromBulkIonDensity2()

    else if(trim(ionic_liquid_name) == "C2MIM+_TFSI-_model2") then
       call SetC2MIN_TFSIBeadDensityFromBulkIonDensity2()

    else if(trim(ionic_liquid_name) == "C6MIM+_TFSI-_model2") then
       call SetC6MIN_TFSIBeadDensityFromBulkIonDensity2()

    else if(trim(ionic_liquid_name) == "C8MIM+_TFSI-_model2") then
       call SetC8MIN_TFSIBeadDensityFromBulkIonDensity2()

    else if(trim(ionic_liquid_name) == "C10MIM+_TFSI-_model2") then
       call SetC10MIN_TFSIBeadDensityFromBulkIonDensity2()

    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then
       call SetPositiveMinusBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralMinusSpheres") then
       call SetPositiveNeutralMinusBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then
       call SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then
       call SetDimerDoubleDimerBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "Pentamers") then
       call SetPentamersBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "Heptamers") then
       call SetHeptamersBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "Heptamer_SingleSphere") then
       call SetHeptamer_SingleSphereBeadDensityFromBulkIonDensity()


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

    positive_oligomer_charge = 0.0_dp
    negative_oligomer_charge = 0.0_dp

  end subroutine SetSingleNeutralSphereBeadDensityFromBulkIonDensity

  subroutine SetSinglePositiveSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = 0.0_dp

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = 0.0_dp

  end subroutine SetSinglePositiveSphereBeadDensityFromBulkIonDensity

  subroutine SetSingleNegativeSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = bulk_density

    positive_oligomer_charge = 0.0_dp
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetSingleNegativeSphereBeadDensityFromBulkIonDensity

  subroutine SetNeutralDimerDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = 2.0_dp * bulk_density
    bulk_density_negative_beads = 0.0_dp

    positive_oligomer_charge = 0.0_dp
    negative_oligomer_charge = 0.0_dp

  end subroutine SetNeutralDimerDensityFromBulkIonDensity

  subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 5.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 5.0_dp * negative_bead_charge

    r_cation = 10
    r_anion = 5
    num_end_monomers_cation = 4
    num_end_monomers_anion = 4

    n_plus_cation_end_bulk = 2.0_dp * bulk_density
    n_neutral_cation_end_bulk = 2.0_dp * bulk_density
    n_minus_cation_end_bulk = 0.0_dp

    n_plus_cation_nonend_bulk = 3.0_dp * bulk_density
    n_neutral_cation_nonend_bulk = 3.0_dp * bulk_density
    n_minus_cation_nonend_bulk = 0.0_dp

    n_plus_anion_end_bulk = 0.0_dp 
    n_neutral_anion_end_bulk = 0.0_dp
    n_minus_anion_end_bulk = 4.0_dp * bulk_density


  end subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetC2MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 3.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 5.0_dp * negative_bead_charge

  end subroutine SetC2MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetC6MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 7.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 5.0_dp * negative_bead_charge

  end subroutine SetC6MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetC8MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 9.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 5.0_dp * negative_bead_charge

  end subroutine SetC8MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetC10MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 11.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 5.0_dp * negative_bead_charge

  end subroutine SetC10MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetC4MIN_TFSIBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 19.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetC4MIN_TFSIBeadDensityFromBulkIonDensity

  subroutine SetC2MIN_TFSIBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 17.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetC2MIN_TFSIBeadDensityFromBulkIonDensity

  subroutine SetC6MIN_TFSIBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 21.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetC6MIN_TFSIBeadDensityFromBulkIonDensity

  subroutine SetC8MIN_TFSIBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 23.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetC8MIN_TFSIBeadDensityFromBulkIonDensity

  subroutine SetC10MIN_TFSIBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 25.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetC10MIN_TFSIBeadDensityFromBulkIonDensity


  subroutine SetC4MIN_TFSIBeadDensityFromBulkIonDensity2()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 16.0_dp * bulk_density
    bulk_density_negative_beads = 4.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 4.0_dp * negative_bead_charge

  end subroutine SetC4MIN_TFSIBeadDensityFromBulkIonDensity2


  subroutine SetC2MIN_TFSIBeadDensityFromBulkIonDensity2()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 14.0_dp * bulk_density
    bulk_density_negative_beads = 4.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 4.0_dp * negative_bead_charge

  end subroutine SetC2MIN_TFSIBeadDensityFromBulkIonDensity2


  subroutine SetC6MIN_TFSIBeadDensityFromBulkIonDensity2()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 18.0_dp * bulk_density
    bulk_density_negative_beads = 4.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 4.0_dp * negative_bead_charge

  end subroutine SetC6MIN_TFSIBeadDensityFromBulkIonDensity2


  subroutine SetC8MIN_TFSIBeadDensityFromBulkIonDensity2()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 20.0_dp * bulk_density
    bulk_density_negative_beads = 4.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 4.0_dp * negative_bead_charge

  end subroutine SetC8MIN_TFSIBeadDensityFromBulkIonDensity2


  subroutine SetC10MIN_TFSIBeadDensityFromBulkIonDensity2()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 22.0_dp * bulk_density
    bulk_density_negative_beads = 4.0_dp * bulk_density

    positive_oligomer_charge = 5.0_dp * positive_bead_charge
    negative_oligomer_charge = 4.0_dp * negative_bead_charge

  end subroutine SetC10MIN_TFSIBeadDensityFromBulkIonDensity2

  subroutine SetPositiveMinusBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = bulk_density

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetPositiveMinusBeadDensityFromBulkIonDensity

  subroutine SetPositiveNeutralMinusBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = bulk_density

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetPositiveNeutralMinusBeadDensityFromBulkIonDensity


  subroutine SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = bulk_density

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity

  subroutine SetDimerDoubleDimerBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 2.0_dp * bulk_density
    bulk_density_neutral_beads = 2.0_dp * bulk_density
    bulk_density_negative_beads = 2.0_dp * bulk_density

    positive_oligomer_charge = 2.0_dp * positive_bead_charge
    negative_oligomer_charge = 2.0_dp * negative_bead_charge

  end subroutine SetDimerDoubleDimerBeadDensityFromBulkIonDensity

  subroutine SetPentamersBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 1.0_dp * bulk_density
    bulk_density_neutral_beads = 8.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = negative_bead_charge
    
  end subroutine SetPentamersBeadDensityFromBulkIonDensity

  subroutine SetHeptamersBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 1.0_dp * bulk_density
    bulk_density_neutral_beads = 12.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetHeptamersBeadDensityFromBulkIonDensity
  
  subroutine SetHeptamer_SingleSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 1.0_dp * bulk_density
    bulk_density_neutral_beads = 6.0_dp * bulk_density
    bulk_density_negative_beads = 1.0_dp * bulk_density

    positive_oligomer_charge = positive_bead_charge
    negative_oligomer_charge = negative_bead_charge

  end subroutine SetHeptamer_SingleSphereBeadDensityFromBulkIonDensity
  

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
