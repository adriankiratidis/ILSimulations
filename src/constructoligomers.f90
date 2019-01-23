!Contains routines to contstruct the various types of oligomers
!that we may wish to calculate.  These routines are present in the module
!as they need to be called in two different places in the code.
!Once during iteration in order to calculate the updated density,
!and then again post density calculation in order to calculate the
!ideal chain contribution to the free energy.
module constructoligomers
  use kinds
  use parameters
  use helpers
  use integratephispherical
  use integratezcylindrical
  use normalisation
  use lambdas
  use charge
  implicit none
  private

  public :: UpdateDensities

  public :: calculate_single_neutral_sphere_ideal_chain_term
  public :: calculate_neutral_dimers_ideal_chain_term
  public :: calculate_PositiveMinusSpheres_ideal_chain_term
  public :: calculate_PositiveNeutralMinusSpheres_ideal_chain_term
  public :: calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term
  public :: calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term
  public :: calculate_C4MIMBF4_ideal_chain_term
  public :: calculate_C2MIMBF4_ideal_chain_term
  public :: calculate_C6MIMBF4_ideal_chain_term
  public :: calculate_C8MIMBF4_ideal_chain_term
  public :: calculate_C10MIMBF4_ideal_chain_term


  public :: calculate_C4MIMTFSI_ideal_chain_term_model1
  public :: calculate_C4MIMTFSI_ideal_chain_term_model2
  
  public :: calculate_C2MIMTFSI_ideal_chain_term_model1
  public :: calculate_C6MIMTFSI_ideal_chain_term_model1
  public :: calculate_C8MIMTFSI_ideal_chain_term_model1
  public :: calculate_C10MIMTFSI_ideal_chain_term_model1

  public :: calculate_C2MIMTFSI_ideal_chain_term_model2
  public :: calculate_C6MIMTFSI_ideal_chain_term_model2
  public :: calculate_C8MIMTFSI_ideal_chain_term_model2
  public :: calculate_C10MIMTFSI_ideal_chain_term_model2
  
  public :: calculate_chem_potential_term_neutral_spheres
  public :: calculate_chem_potential_term_neutral_dimers
  public :: calculate_chem_potential_C4MIMBF4

  !public :: calculate_chem_potential_C2MIMBF4
  !public :: calculate_chem_potential_C6MIMBF4
  !public :: calculate_chem_potential_C8MIMBF4
  !public :: calculate_chem_potential_C10MIMBF4

  public :: calculate_chem_potential_C4MIMTFSI_model1
  public :: calculate_chem_potential_C4MIMTFSI_model2
  !public :: calculate_chem_potential_C2MIMTFSI_model1
  !public :: calculate_chem_potential_C6MIMTFSI_model1
  !public :: calculate_chem_potential_C8MIMTFSI_model1
  !public :: calculate_chem_potential_C10MIMTFSI_model1
  public :: calculate_chem_potential_PositiveMinusSpheres
  public :: calculate_chem_potential_PositiveNeutralMinusSpheres
  public :: calculate_chem_potential_PositiveNeutralDimerMinusSpheres
  public :: calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer

  !Private subroutines
  !UpdateSinglePositiveSphereDensity
  !UpdateSingleNeutralSphereDensity
  !UpdateSinglePositiveSphereDensity
  !UpdateSinglePositiveNeutralMinusSphereDensities
  !UpdateC4MIMBF4PositiveBeadDensities
  !UpdateC2MIMPositiveBeadDensities
  !UpdateC6MIMPositiveBeadDensities
  !UpdateC8MIMPositiveBeadDensities
  !UpdateC10MIMPositiveBeadDensities
  !UpdateC4MIMBF4NeutralBeadDensities
  !UpdateC4MIMBF4NegativeBeadDensities
  !UpdateC4MIMTFSINeutralBeadDensities_model1
  !UpdateC4MIMTFSINegativeBeadDensities_model1
  !UpdateC4MIMTFSINeutralBeadDensities_model2
  !UpdateC4MIMTFSINegativeBeadDensities_model2


contains

  subroutine UpdateDensities(n1, n2, n3, lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, Donnan_potential, iteration, abort_now)
    real(dp), dimension(:), intent(in) :: n1
    real(dp), dimension(:), intent(in) :: n2
    real(dp), dimension(:), intent(in) :: n3

    real(dp), dimension(:), intent(in) :: lambda1
    real(dp), dimension(:), intent(out) :: n1_updated

    real(dp), dimension(:), intent(in), optional :: lambda2
    real(dp), dimension(:), intent(out), optional :: n2_updated

    real(dp), dimension(:), intent(in), optional :: lambda3
    real(dp), dimension(:), intent(out), optional :: n3_updated

    real(dp), dimension(:), intent(in) :: lambda_hs_end
    real(dp), dimension(:), intent(out) :: n_hs_end

    real(dp), dimension(:), intent(in) :: lambda_hs_nonend
    real(dp), dimension(:), intent(out) :: n_hs_nonend
    
    real(dp) :: Donnan_potential
    integer :: iteration
    logical :: abort_now

    real(dp) :: intermediate
    integer :: ij

    real(dp) :: Donnan_potential_previous

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNeutralSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSinglePositiveSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNegativeSphereDensity(lambda1, n1_updated)


    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveMinusSphereDensities(lambda1, n1_updated, lambda3, n3_updated)
       end if

    else if(trim(ionic_liquid_name) == "PositiveNeutralMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateSinglePositiveNeutralMinusSphereDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)


          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous



       end if



       ! else if(trim(ionic_liquid_name) == "PositiveNeutralMinusSpheres") then

       !    if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
       !         (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
       !       print *, "constructoligomers.90: UpdateDensities:"
       !       print *, "Trying to Update positive, neutral and negative spheres"
       !       print *, "but haven't been passed the appropriate number of arguments"
       !       print *, "coding error...aborting..."
       !       call abort()
       !    else
       !       call UpdatePositiveNeutralMinusSphereDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
       !    end if


    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveNeutralDimerMinusSphereDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
       end if

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveNeutralDoubleDimerMinusDimerDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous


       end if

    else if(trim(ionic_liquid_name) == "NeutralDimers") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing Neutral Dimers should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       else
          !call UpdateSingleNeutralSphereDensity(lambda1, n1_updated)
          call UpdateNeutralDimerDensity(lambda1, n1_updated)
       end if

    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          
          call UpdateC4MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
          call UpdateC4MIMBF4NeutralBeadDensities(lambda1, lambda2, n2_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
          call UpdateC4MIMBF4NegativeBeadDensities(lambda3, n3_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous


       end if


    else if(trim(ionic_liquid_name) == "C2MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC2MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC2MIMBF4NeutralBeadDensities(lambda1, lambda2, n2_updated)
          call UpdateC4MIMBF4NegativeBeadDensities(lambda3, n3_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous


       end if

    else if(trim(ionic_liquid_name) == "C6MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC6MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC6MIMBF4NeutralBeadDensities(lambda1, lambda2, n2_updated)
          call UpdateC4MIMBF4NegativeBeadDensities(lambda3, n3_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous


       end if

    else if(trim(ionic_liquid_name) == "C8MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC8MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC8MIMBF4NeutralBeadDensities(lambda1, lambda2, n2_updated)
          call UpdateC4MIMBF4NegativeBeadDensities(lambda3, n3_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous


       end if

    else if(trim(ionic_liquid_name) == "C10MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC10MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC10MIMBF4NeutralBeadDensities(lambda1, lambda2, n2_updated)
          call UpdateC4MIMBF4NegativeBeadDensities(lambda3, n3_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n2_updated = n2_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous


       end if



    else if(trim(ionic_liquid_name) == "C2MIM+_TFSI-_model1") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC2MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model1(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC2MIMTFSINeutralBeadDensities_model1(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

    else if(trim(ionic_liquid_name) == "C6MIM+_TFSI-_model1") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC6MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model1(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC6MIMTFSINeutralBeadDensities_model1(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

    else if(trim(ionic_liquid_name) == "C8MIM+_TFSI-_model1") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC8MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model1(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC8MIMTFSINeutralBeadDensities_model1(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

    else if(trim(ionic_liquid_name) == "C10MIM+_TFSI-_model1") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC10MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model1(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC10MIMTFSINeutralBeadDensities_model1(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if


    else if(trim(ionic_liquid_name) == "C4MIM+_TFSI-_model1") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC4MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
          call UpdateC4MIMTFSINegativeBeadDensities_model1(lambda2, lambda3, n3_updated)


          ! if(iteration < 1000) then

          !    do ij = 1, size(n1)
          !       if(n1_updated(ij) .gt. n1(ij)/exp(beta*Donnan_potential*positive_oligomer_charge)) then
          !          intermediate = 2*n1(ij)/exp(beta*Donnan_potential*positive_oligomer_charge) &
          !               - ((n1(ij)/exp(beta*Donnan_potential*positive_oligomer_charge))**2)/n1_updated(ij)
          !          if(intermediate .lt. n1_updated(ij)) then
          !             n1_updated(ij) = intermediate
          !          end if
          !       end if
          !    end do

          !    do ij = 1, size(n3)
          !       if(n3_updated(ij) .gt. n3(ij)/exp(beta*Donnan_potential*negative_oligomer_charge)) then
          !          intermediate = 2*n3(ij)/exp(beta*Donnan_potential*negative_oligomer_charge) &
          !               - ((n3(ij)/exp(beta*Donnan_potential*negative_oligomer_charge))**2)/n3_updated(ij)
          !          if(intermediate .lt. n3_updated(ij)) then
          !             n3_updated(ij) = intermediate
          !          end if
          !       end if
          !    end do

          ! end if

          !if(iteration > 1) then
          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)
          !end if



          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)
          !Donnan_potential = 0.0_dp


          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous

          call UpdateC4MIMTFSINeutralBeadDensities_model1(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)


       end if

    else if(trim(ionic_liquid_name) == "C4MIM+_TFSI-_model2") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC4MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
          call UpdateC4MIMTFSINegativeBeadDensities_model2(lambda2, lambda3, n3_updated)

          !if(iteration > 1) then
          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)
          !end if

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)


          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous

          call UpdateC4MIMTFSINeutralBeadDensities_model2(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

       
    else if(trim(ionic_liquid_name) == "C2MIM+_TFSI-_model2") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC2MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model2(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC2MIMTFSINeutralBeadDensities_model2(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

    else if(trim(ionic_liquid_name) == "C6MIM+_TFSI-_model2") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC6MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model2(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC6MIMTFSINeutralBeadDensities_model2(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

    else if(trim(ionic_liquid_name) == "C8MIM+_TFSI-_model2") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC8MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model2(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC8MIMTFSINeutralBeadDensities_model2(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

    else if(trim(ionic_liquid_name) == "C10MIM+_TFSI-_model2") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC10MIMBF4PositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMTFSINegativeBeadDensities_model2(lambda2, lambda3, n3_updated)

          Donnan_potential_previous = Donnan_potential
          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          call CalculateDonnanPotential(n1_updated, n3_updated, Donnan_potential, abort_now)

          n1_updated = n1_updated*exp(beta*Donnan_potential*positive_oligomer_charge)
          n3_updated = n3_updated*exp(beta*Donnan_potential*negative_oligomer_charge)

          Donnan_potential = Donnan_potential + Donnan_potential_previous
          call UpdateC10MIMTFSINeutralBeadDensities_model2(lambda1, lambda2, lambda3, n2_updated, Donnan_potential)

       end if

       
    else

       print *, "constructoligomers.f90: UpdateDensities: "
       print *, "Unsupported ionic_liquid_name value of ", trim(ionic_liquid_name)
       print *, "...aborting..."
       call abort()
    end if

    !Imposes charge neutrality, by solving for the Donnan potential.
    if(present(n3_updated) .and. trim(ionic_liquid_name) /= "C4MIM+_TFSI-_model1" &
         .and. trim(ionic_liquid_name) /= "C4MIM+_TFSI-_model2" .and. trim(ionic_liquid_name) /= "C4MIM_BF4-" &
         .and. trim(ionic_liquid_name) /= "C2MIM+_TFSI-_model1" .and. trim(ionic_liquid_name) /= "C6MIM+_TFSI-_model1" &
         .and. trim(ionic_liquid_name) /= "C8MIM+_TFSI-_model1" .and. trim(ionic_liquid_name) /= "C10MIM+_TFSI-_model1" &
         .and. trim(ionic_liquid_name) /= "PositiveNeutralDoubleDimerMinusDimer" .and. trim(ionic_liquid_name) /= "C2MIM_BF4-" &
         .and. trim(ionic_liquid_name) /= "C6MIM_BF4-" .and. trim(ionic_liquid_name) /= "C8MIM_BF4-" &
         .and. trim(ionic_liquid_name) /= "C10MIM_BF4-" .and. trim(ionic_liquid_name) /= "C2MIM+_TFSI-_model2" &
         .and. trim(ionic_liquid_name) /= "C6MIM+_TFSI-_model2" .and. trim(ionic_liquid_name) /= "C8MIM+_TFSI-_model2" &
         .and. trim(ionic_liquid_name) /= "C10MIM+_TFSI-_model2" &
         ) then !we have positive, neutral and negative spheres all present
       !Note: if only doing a single charge species then this conditional call needs updating.
       !call ImposeChargeNeutrality(n1_updated, n2_updated, n3_updated, Donnan_potential, abort_now)
    end if

  end subroutine UpdateDensities

  !****************************************************************
  !*****START OF UPDATE DENSITY ROUTINES
  !****************************************************************

  subroutine UpdateSinglePositiveSphereDensity(lambda_plus, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated

    if(size(lambda_plus) /= size(n_plus_updated)) then
       print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
       print *, "Size mismatch.  size(lambda_plus) /= size(n_plus_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_plus_updated = bulk_density_positive_beads * exp(lambda_plus)
    end if

    call setNonCalculatedRegionToZero(n_plus_updated)
    !call RenormaliseToBulkDensity(n_plus_updated, "n+")

  end subroutine UpdateSinglePositiveSphereDensity


  subroutine UpdateSingleNeutralSphereDensity(lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    if(size(lambda_neutral) /= size(n_neutral_updated)) then
       print *, "constructoligomers.f90: UpdateSingleNeutralSphereDensity: "
       print *, "Size mismatch.  size(lambda_neutral) /= size(n_plus_neutral)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_neutral_updated = bulk_density_neutral_beads * exp(lambda_neutral)
    end if

    call setNonCalculatedRegionToZero(n_neutral_updated)
    !call RenormaliseToBulkDensity(n_neutral_updated, "n0")

  end subroutine UpdateSingleNeutralSphereDensity


  subroutine UpdateSingleNegativeSphereDensity(lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    if(size(lambda_minus) /= size(n_minus_updated)) then
       print *, "constructoligomers.f90: UpdateSingleNegativeSphereDensity: "
       print *, "Size mismatch.  size(lambda_plus) /= size(n_plus_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_minus_updated = bulk_density_negative_beads * exp(lambda_minus)
    end if

    call setNonCalculatedRegionToZero(n_minus_updated)
    !call RenormaliseToBulkDensity(n_minus_updated, "n-")

  end subroutine UpdateSingleNegativeSphereDensity


  subroutine UpdateSinglePositiveNeutralMinusSphereDensities(lambda_plus, n_plus_updated, &
       lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: UpdateSinglePositiveNeutralMinusSphereDensities:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    call UpdateSinglePositiveSphereDensity(lambda_plus, n_plus_updated)
    call UpdateSingleNeutralSphereDensity(lambda_neutral, n_neutral_updated)
    call UpdateSingleNegativeSphereDensity(lambda_minus, n_minus_updated)

  end subroutine UpdateSinglePositiveNeutralMinusSphereDensities


  subroutine UpdateNeutralDimerDensity(lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), allocatable :: c1

    ! First check the input variables are the same size
    if( (size(lambda_neutral) == size(n_neutral_updated)) ) then
       allocate(c1(size(lambda_neutral)))
    else
       print *, "constructoligomers.f90: UpdateNeutralDimerDensity:"
       print *, "Size mismatch. size(lambda_neutral) /= size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    !c1 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral))
    c1 = integrate_phi_spherical(exp(lambda_neutral))
    !c1 = 0.5_dp * (exp(lambda_neutral))

    !Note the factor of 2 is present as this formulation calculates the density of
    !an individual bead.  By symmetry we need both the beads are the same and to
    !get the total bead density we must add them.  Hence the factor of 2.
    n_neutral_updated = 2.0_dp * bulk_density * exp(lambda_neutral) * c1

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c1)

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateNeutralDimerDensity


  subroutine UpdatePositiveMinusSphereDensities(lambda_plus, n_plus_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    integer :: ij

    !print *, "lambda_plus = ", lambda_plus
    !print *, "lambda_minus = ", lambda_minus

    n_plus_updated = bulk_density * exp(lambda_plus)
    n_minus_updated = bulk_density * exp(lambda_minus)

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
       n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    end do

    !print * , "n_plus_updated = ", n_plus_updated
    !print * , "n_minus_updated = ", n_minus_updated
    !call abort()
    call setNonCalculatedRegionToZero(n_plus_updated)
    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdatePositiveMinusSphereDensities


  subroutine UpdatePositiveNeutralDimerMinusSphereDensities(lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), allocatable :: c1, c2
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated)) ) then
       allocate(c1(size(lambda_plus)))
       allocate(c2(size(lambda_plus)))
    else
       print *, "constructoligomers.f90: UpdatePositiveNeutralDimerMinusSphereDensities:"
       print *, "Size mismatch. The following expression is false."
       print *, "(size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated))"
       print *, "...aborting..."
       call abort()
    end if

    c1 = integrate_phi_spherical(exp(lambda_plus))
    c2 = integrate_phi_spherical(exp(lambda_neutral))

    n_plus_updated = bulk_density * exp(lambda_plus) * c2
    n_neutral_updated = bulk_density * exp(lambda_neutral) * c1

    n_minus_updated = bulk_density * exp(lambda_minus)

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
       n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    end do

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c1, c2)

    call setNonCalculatedRegionToZero(n_plus_updated)
    call setNonCalculatedRegionToZero(n_neutral_updated)
    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdatePositiveNeutralDimerMinusSphereDensities


  subroutine UpdatePositiveNeutralDoubleDimerMinusDimerDensities(lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), allocatable :: c1, c2, c3, c4
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated)) ) then
       allocate(c1(size(lambda_plus)))
       allocate(c2(size(lambda_plus)))
       allocate(c3(size(lambda_plus)))
       allocate(c4(size(lambda_plus)))       
    else
       print *, "constructoligomers.f90: UpdatePositiveNeutralDimerMinusSphereDensities:"
       print *, "Size mismatch. The following expression is false."
       print *, "(size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated))"
       print *, "...aborting..."
       call abort()
    end if

    c1 = integrate_phi_spherical(exp(lambda_minus))
    n_minus_updated = 2.0_dp * bulk_density * exp(lambda_minus) * c1

    c4 = integrate_phi_spherical(exp(lambda_neutral))
    c3 = integrate_phi_spherical(exp(lambda_neutral) * c4)
    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c1 = integrate_phi_spherical(exp(lambda_plus))

    n_plus_updated = bulk_density * ((exp(lambda_plus) * c2) + (exp(lambda_plus) * c3 * c1))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c1)
    c3 = integrate_phi_spherical(exp(lambda_neutral) * c2)

    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c3) + (exp(lambda_neutral) * c2 * c4) )



    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
    end do

    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    do ij = 1, (size(n_minus_updated) - 1)/2
       n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    end do




    deallocate(c1, c2, c3, c4)

    call setNonCalculatedRegionToZero(n_plus_updated)
    call setNonCalculatedRegionToZero(n_neutral_updated)
    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdatePositiveNeutralDoubleDimerMinusDimerDensities


  subroutine UpdateC4MIMBF4PositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), intent(in) :: lambda_hs_end
    real(dp), dimension(:), intent(out) :: n_hs_end
    real(dp), dimension(:), intent(in) :: lambda_hs_nonend
    real(dp), dimension(:), intent(out) :: n_hs_nonend

    real(dp), dimension(size(lambda_plus)) :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_plus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = integrate_phi_spherical(exp(lambda_plus + lambda_hs_end))

    c8c1 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_end))

    c7 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c8c1)

    c6 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c5)

    c3p = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c4 * c9c10 * c9c10)

    c2 = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c8c1)

    c3pp = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c4 * c2 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus + lambda_hs_nonend) * c8c1 * c3p) + (exp(lambda_plus + lambda_hs_nonend) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus + lambda_hs_nonend) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus + lambda_hs_end) * c3pp)) )

    n_hs_end = n_hs_end +  bulk_density * ( (2.0_dp * (exp(lambda_plus + lambda_hs_end) * c3pp)) )
    
    n_hs_nonend = n_hs_nonend + bulk_density * ( (exp(lambda_plus + lambda_hs_nonend) * c8c1 * c3p) + (exp(lambda_plus + lambda_hs_nonend) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus + lambda_hs_nonend) * c3ppp * c5) )



    !n_plus_updated = 0.0_dp

    !do ij = 1, (size(n_plus_updated) - 1)/2
    !   n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
    !end do

    !call RenormaliseToBulkDensity(n_plus_updated, "n+")
    call setNonCalculatedRegionToZero(n_plus_updated)

  end subroutine UpdateC4MIMBF4PositiveBeadDensities


  subroutine UpdateC2MIMBF4PositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(size(lambda_plus)) :: c6c1, c9c10, c5, c4, c2, c3p, c3pp, c3ppp
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_plus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c6c1 = integrate_phi_spherical(exp(lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6c1)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c6c1)

    c3pp = integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus) * c6c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    !n_plus_updated = 0.0_dp

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
    end do

    !call RenormaliseToBulkDensity(n_plus_updated, "n+")
    call setNonCalculatedRegionToZero(n_plus_updated)

  end subroutine UpdateC2MIMBF4PositiveBeadDensities


  subroutine UpdateC6MIMBF4PositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(size(lambda_plus)) :: c12c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp, c11, c8
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_plus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c12c1 = integrate_phi_spherical(exp(lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12c1)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c12c1)

    c3pp = integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus) * c12c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    !n_plus_updated = 0.0_dp

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
    end do

    !call RenormaliseToBulkDensity(n_plus_updated, "n+")
    call setNonCalculatedRegionToZero(n_plus_updated)

  end subroutine UpdateC6MIMBF4PositiveBeadDensities

  subroutine UpdateC8MIMBF4PositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(size(lambda_plus)) :: c14c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp, c11, c8, c12, c13
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_plus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c14c1 = integrate_phi_spherical(exp(lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14c1)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c14c1)

    c3pp = integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus) * c14c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    !n_plus_updated = 0.0_dp

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
    end do

    !call RenormaliseToBulkDensity(n_plus_updated, "n+")
    call setNonCalculatedRegionToZero(n_plus_updated)

  end subroutine UpdateC8MIMBF4PositiveBeadDensities


  subroutine UpdateC10MIMBF4PositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(size(lambda_plus)) :: c16c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp, c11, c8, c12, c13, c14, c15
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_plus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c16c1 = integrate_phi_spherical(exp(lambda_neutral))

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16c1)

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c16c1)

    c3pp = integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus) * c16c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    !n_plus_updated = 0.0_dp

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
    end do

    !call RenormaliseToBulkDensity(n_plus_updated, "n+")
    call setNonCalculatedRegionToZero(n_plus_updated)

  end subroutine UpdateC10MIMBF4PositiveBeadDensities


  subroutine UpdateC4MIMBF4NeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), intent(in) :: lambda_hs_end
    real(dp), dimension(:), intent(out) :: n_hs_end
    real(dp), dimension(:), intent(in) :: lambda_hs_nonend
    real(dp), dimension(:), intent(out) :: n_hs_nonend

    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2 !contributions also used in +ve beads.
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p !extra contributions for -ve beads.

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .and. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    c9c10 = integrate_phi_spherical(exp(lambda_plus + lambda_hs_end))

    c8c1 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_end))

    c7 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c8c1)

    c6 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c8c1)

    c3p = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus + lambda_hs_nonend) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral + lambda_hs_nonend) * c6p)

    !Calculate the resulting neutral bead densities.
    !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral + lambda_hs_end) * c2p) + (exp(lambda_neutral + lambda_hs_nonend) * c4p * c6) + &
         (exp(lambda_neutral + lambda_hs_nonend) * c5p * c7) + (exp(lambda_neutral + lambda_hs_nonend) * c6p * c8c1) + (exp(lambda_neutral + lambda_hs_end) * c7p) )

    n_hs_end = n_hs_end + bulk_density * ( (exp(lambda_neutral + lambda_hs_end) * c2p) + (exp(lambda_neutral + lambda_hs_end) * c7p) )
    
    n_hs_nonend = n_hs_nonend + bulk_density * ( (exp(lambda_neutral + lambda_hs_nonend) * c4p * c6) + &
         (exp(lambda_neutral + lambda_hs_nonend) * c5p * c7) + (exp(lambda_neutral + lambda_hs_nonend) * c6p * c8c1) )


    ! do ij = 1, (size(n_neutral_updated) - 1)/2
    !    n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    ! end do


    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC4MIMBF4NeutralBeadDensities


  subroutine UpdateC4MIMBF4NegativeBeadDensities(lambda_minus, n_minus_updated, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), intent(in) :: lambda_hs_end
    real(dp), dimension(:), intent(out) :: n_hs_end
    real(dp), dimension(:), intent(in) :: lambda_hs_nonend
    real(dp), dimension(:), intent(out) :: n_hs_nonend

    real(dp), dimension(size(lambda_minus)) :: a1a2a3a4, a5p
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if(size(lambda_minus) /= size(n_minus_updated)) then
       print *, "constructoligomers.f90: UpdateBF4NegativeBeads:"
       print *, "Size mismatch. size(lambda_minus) /= size(n_minus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    
    !Calculate the required contributions for the anion
    a1a2a3a4 = integrate_phi_spherical(exp(lambda_minus + lambda_hs_end))

    
    a5p = integrate_phi_spherical(exp(lambda_minus + lambda_hs_nonend) * (a1a2a3a4 ** 3.0_dp))

    
    !Calculate the resulting negative bead densities.
    !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
    n_minus_updated = bulk_density * ( 4.0_dp*(exp(lambda_minus + lambda_hs_end) * a5p) + (exp(lambda_minus + lambda_hs_nonend) * (a1a2a3a4**4.0_dp)) )

    n_hs_end = n_hs_end + bulk_density * ( 4.0_dp*(exp(lambda_minus + lambda_hs_end) * a5p) )
    n_hs_nonend = n_hs_nonend + bulk_density * ( (exp(lambda_minus + lambda_hs_nonend) * (a1a2a3a4**4.0_dp)) )

    
    ! do ij = 1, (size(n_minus_updated) - 1)/2
    !    n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    ! end do

    call setNonCalculatedRegionToZero(n_minus_updated)


  end subroutine UpdateC4MIMBF4NegativeBeadDensities


  subroutine UpdateC4MIMTFSINeutralBeadDensities_model1(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c8c1 = integrate_phi_spherical(exp(lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8c1)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c8c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1) + (exp(lambda_neutral) * c7p) )


    !Now add the density contribution to n_neutral from the anion.
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = F

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3) 

    NCOOCF3 = integrate_phi_spherical(exp(lambda_minus) * COOCF3) 

    F3CCONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3) 

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (4.0_dp * bulk_density * exp(lambda_neutral) * F3CCONCOOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now explicitly ensure symmetry
    !do ij = 1, (size(n_neutral_updated) - 1)/2
    !   n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    !end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC4MIMTFSINeutralBeadDensities_model1

  subroutine UpdateC2MIMTFSINeutralBeadDensities_model1(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c6c1, c7, c6, c5, c4, c2
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c6c1 = integrate_phi_spherical(exp(lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6c1)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c6c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)


    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6c1) + &
         (exp(lambda_neutral) * c5p) )


    !Now add the density contribution to n_neutral from the anion.
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = F

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3) 

    NCOOCF3 = integrate_phi_spherical(exp(lambda_minus) * COOCF3) 

    F3CCONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3) 

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (4.0_dp * bulk_density * exp(lambda_neutral) * F3CCONCOOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC2MIMTFSINeutralBeadDensities_model1


  subroutine UpdateC6MIMTFSINeutralBeadDensities_model1(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c12c1, c7, c6, c5, c4, c2, c8, c11
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c12c1 = integrate_phi_spherical(exp(lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12c1)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c12c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12c1) + &
         (exp(lambda_neutral) * c11p) )


    !Now add the density contribution to n_neutral from the anion.
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = F

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3) 

    NCOOCF3 = integrate_phi_spherical(exp(lambda_minus) * COOCF3) 

    F3CCONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3) 

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (4.0_dp * bulk_density * exp(lambda_neutral) * F3CCONCOOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC6MIMTFSINeutralBeadDensities_model1


  subroutine UpdateC8MIMTFSINeutralBeadDensities_model1(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c14c1, c7, c6, c5, c4, c2, c8, c11, c12, c13
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p, c12p, c13p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c14c1 = integrate_phi_spherical(exp(lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14c1)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c14c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    c12p = integrate_phi_spherical(exp(lambda_neutral) * c11p)

    c13p = integrate_phi_spherical(exp(lambda_neutral) * c12p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12) + &
         (exp(lambda_neutral) * c11p * c13) + (exp(lambda_neutral) * c12p * c14c1) + (exp(lambda_neutral) * c13p) )


    !Now add the density contribution to n_neutral from the anion.
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = F

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3) 

    NCOOCF3 = integrate_phi_spherical(exp(lambda_minus) * COOCF3) 

    F3CCONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3) 

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (4.0_dp * bulk_density * exp(lambda_neutral) * F3CCONCOOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC8MIMTFSINeutralBeadDensities_model1


  subroutine UpdateC10MIMTFSINeutralBeadDensities_model1(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c16c1, c7, c6, c5, c4, c2, c8, c11, c12, c13, c14, c15
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p, c12p, c13p, c14p, c15p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c16c1 = integrate_phi_spherical(exp(lambda_neutral))

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16c1)

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c16c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    c12p = integrate_phi_spherical(exp(lambda_neutral) * c11p)

    c13p = integrate_phi_spherical(exp(lambda_neutral) * c12p)

    c14p = integrate_phi_spherical(exp(lambda_neutral) * c13p)

    c15p = integrate_phi_spherical(exp(lambda_neutral) * c14p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12) + &
         (exp(lambda_neutral) * c11p * c13) + (exp(lambda_neutral) * c12p * c14) + (exp(lambda_neutral) * c13p * c15) + (exp(lambda_neutral) * c14p * c16c1) + &
         (exp(lambda_neutral) * c15p) )

    !Now add the density contribution to n_neutral from the anion.
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = F

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3) 

    NCOOCF3 = integrate_phi_spherical(exp(lambda_minus) * COOCF3) 

    F3CCONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3) 

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (4.0_dp * bulk_density * exp(lambda_neutral) * F3CCONCOOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC10MIMTFSINeutralBeadDensities_model1

  
  subroutine UpdateC2MIMTFSINeutralBeadDensities_model2(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c6c1, c7, c6, c5, c4, c2
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p
    
    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O
    
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c6c1 = integrate_phi_spherical(exp(lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6c1)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c6c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)


    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6c1) + &
         (exp(lambda_neutral) * c5p) )

    
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = integrate_phi_spherical(exp(lambda_minus))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3)

    NCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * COOCF3)

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (1.0_dp * bulk_density * exp(lambda_neutral) * COOCF3 * COOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))
    
    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC2MIMTFSINeutralBeadDensities_model2


  subroutine UpdateC6MIMTFSINeutralBeadDensities_model2(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c12c1, c7, c6, c5, c4, c2, c8, c11
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p
    
    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O
    
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c12c1 = integrate_phi_spherical(exp(lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12c1)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c12c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12c1) + &
         (exp(lambda_neutral) * c11p) )


    F = integrate_phi_spherical(exp(lambda_neutral))
    O = integrate_phi_spherical(exp(lambda_minus))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3)

    NCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * COOCF3)

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (1.0_dp * bulk_density * exp(lambda_neutral) * COOCF3 * COOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))


    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC6MIMTFSINeutralBeadDensities_model2


  subroutine UpdateC8MIMTFSINeutralBeadDensities_model2(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c14c1, c7, c6, c5, c4, c2, c8, c11, c12, c13
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p, c12p, c13p
    
    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c14c1 = integrate_phi_spherical(exp(lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14c1)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c14c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    c12p = integrate_phi_spherical(exp(lambda_neutral) * c11p)

    c13p = integrate_phi_spherical(exp(lambda_neutral) * c12p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12) + &
         (exp(lambda_neutral) * c11p * c13) + (exp(lambda_neutral) * c12p * c14c1) + (exp(lambda_neutral) * c13p) )

    
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = integrate_phi_spherical(exp(lambda_minus))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3)

    NCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * COOCF3)

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (1.0_dp * bulk_density * exp(lambda_neutral) * COOCF3 * COOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    
    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC8MIMTFSINeutralBeadDensities_model2


  subroutine UpdateC10MIMTFSINeutralBeadDensities_model2(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c16c1, c7, c6, c5, c4, c2, c8, c11, c12, c13, c14, c15
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p, c12p, c13p, c14p, c15p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c16c1 = integrate_phi_spherical(exp(lambda_neutral))

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16c1)

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c16c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    c12p = integrate_phi_spherical(exp(lambda_neutral) * c11p)

    c13p = integrate_phi_spherical(exp(lambda_neutral) * c12p)

    c14p = integrate_phi_spherical(exp(lambda_neutral) * c13p)

    c15p = integrate_phi_spherical(exp(lambda_neutral) * c14p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12) + &
         (exp(lambda_neutral) * c11p * c13) + (exp(lambda_neutral) * c12p * c14) + (exp(lambda_neutral) * c13p * c15) + (exp(lambda_neutral) * c14p * c16c1) + &
         (exp(lambda_neutral) * c15p) )

    
    F = integrate_phi_spherical(exp(lambda_neutral))
    O = integrate_phi_spherical(exp(lambda_minus))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3)

    NCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * COOCF3)

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (1.0_dp * bulk_density * exp(lambda_neutral) * COOCF3 * COOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))


    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC10MIMTFSINeutralBeadDensities_model2

  
  subroutine UpdateC2MIMBF4NeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c6c1, c7, c6, c5, c4, c2
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c6c1 = integrate_phi_spherical(exp(lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6c1)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c6c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)


    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6c1) + &
         (exp(lambda_neutral) * c5p) )


    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC2MIMBF4NeutralBeadDensities


  subroutine UpdateC6MIMBF4NeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    
    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c12c1, c7, c6, c5, c4, c2, c8, c11
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c12c1 = integrate_phi_spherical(exp(lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12c1)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c12c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12c1) + &
         (exp(lambda_neutral) * c11p) )

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC6MIMBF4NeutralBeadDensities


  subroutine UpdateC8MIMBF4NeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c14c1, c7, c6, c5, c4, c2, c8, c11, c12, c13
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p, c12p, c13p

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c14c1 = integrate_phi_spherical(exp(lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14c1)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c14c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    c12p = integrate_phi_spherical(exp(lambda_neutral) * c11p)

    c13p = integrate_phi_spherical(exp(lambda_neutral) * c12p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12) + &
         (exp(lambda_neutral) * c11p * c13) + (exp(lambda_neutral) * c12p * c14c1) + (exp(lambda_neutral) * c13p) )


    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC8MIMBF4NeutralBeadDensities


  subroutine UpdateC10MIMBF4NeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c16c1, c7, c6, c5, c4, c2, c8, c11, c12, c13, c14, c15
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p, c8p, c11p, c12p, c13p, c14p, c15p

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c16c1 = integrate_phi_spherical(exp(lambda_neutral))

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16c1)

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c16c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    c8p = integrate_phi_spherical(exp(lambda_neutral) * c7p)

    c11p = integrate_phi_spherical(exp(lambda_neutral) * c8p)

    c12p = integrate_phi_spherical(exp(lambda_neutral) * c11p)

    c13p = integrate_phi_spherical(exp(lambda_neutral) * c12p)

    c14p = integrate_phi_spherical(exp(lambda_neutral) * c13p)

    c15p = integrate_phi_spherical(exp(lambda_neutral) * c14p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8) + (exp(lambda_neutral) * c7p * c11) + (exp(lambda_neutral) * c8p * c12) + &
         (exp(lambda_neutral) * c11p * c13) + (exp(lambda_neutral) * c12p * c14) + (exp(lambda_neutral) * c13p * c15) + (exp(lambda_neutral) * c14p * c16c1) + &
         (exp(lambda_neutral) * c15p) )

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC10MIMBF4NeutralBeadDensities


  subroutine UpdateC4MIMTFSINegativeBeadDensities_model1(lambda_neutral, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(size(lambda_neutral)) :: F, CF3, CF3OO
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_neutral) /= size(lambda_minus)) .or. (size(lambda_minus) /= size(n_minus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_neutral) /= size(lambda_minus)"
       print *, "or size(lambda_minus) == size(n_minus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    F = integrate_phi_spherical(exp(lambda_neutral))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    CF3OO = integrate_phi_spherical(exp(lambda_neutral) * CF3 * F * F) !Note: F = "O" as both the F and O atoms and modelled as neutral

    n_minus_updated = bulk_density * exp(lambda_minus) * CF3OO * CF3OO

    !Explicitly ensure symmetry
    !do ij = 1, (size(n_minus_updated) - 1)/2
    !   n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    !end do

    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdateC4MIMTFSINegativeBeadDensities_model1

  subroutine UpdateC4MIMTFSINeutralBeadDensities_model2(lambda_plus, lambda_neutral, lambda_minus, n_neutral_updated, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), intent(in) :: Donnan_potential

    !In order to calculate the contribution from the cation.
    real(dp), dimension(size(lambda_plus)) :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2
    real(dp), dimension(size(lambda_plus)) :: c2p, c4p, c5p, c6p, c7p

    !In order to calculate the contribution from the anion.
    real(dp), dimension(size(lambda_plus)) :: F2CCOONCOOCF3, COONCOOCF3, F3CCONCOOCF3, NCOOCF3, COOCF3, CF3, F, O

    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(n_neutral_updated))) then
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !First calculate the contribution to n_neutral from the cation
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c8c1 = integrate_phi_spherical(exp(lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8c1)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c8c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    !Calculate the resulting neutral bead densities from the cation.
    !n_neutral_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * exp(beta*Donnan_potential*positive_oligomer_charge) * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1) + (exp(lambda_neutral) * c7p) )


    F = integrate_phi_spherical(exp(lambda_neutral))
    O = integrate_phi_spherical(exp(lambda_minus))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3)

    NCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * COOCF3)

    COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
         (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
         (1.0_dp * bulk_density * exp(lambda_neutral) * COOCF3 * COOCF3) + &
         (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
         (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now add the density contribution to n_neutral from the anion.

    ! F = integrate_phi_spherical(exp(lambda_neutral))
    ! O = F

    ! CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    ! COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3) 

    ! NCOOCF3 = integrate_phi_spherical(exp(lambda_minus) * COOCF3) 

    ! F3CCONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3) 

    ! COONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * NCOOCF3) 

    ! F2CCOONCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * COONCOOCF3) 

    ! n_neutral_updated = n_neutral_updated + (exp(beta*Donnan_potential*negative_oligomer_charge) *(&
    !      (2.0_dp * bulk_density * exp(lambda_neutral) * NCOOCF3 * O * O * CF3) + &
    !      (4.0_dp * bulk_density * exp(lambda_neutral) * F3CCONCOOCF3) + &
    !      (2.0_dp * bulk_density * exp(lambda_neutral) * COONCOOCF3 * F * F * F) + &
    !      (6.0_dp * bulk_density * exp(lambda_neutral) * F2CCOONCOOCF3)))

    !Now explicitly ensure symmetry
    do ij = 1, (size(n_neutral_updated) - 1)/2
       n_neutral_updated(ij) = n_neutral_updated(size(n_neutral_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateC4MIMTFSINeutralBeadDensities_model2


  subroutine UpdateC4MIMTFSINegativeBeadDensities_model2(lambda_neutral, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(size(lambda_neutral)) :: F, CF3, CF3OO, O, NCOOCF3, NCOOCF3COCF3, COOCF3
    integer :: array_size
    integer :: ij

    ! First check the input variables are the same size
    if( (size(lambda_neutral) /= size(lambda_minus)) .or. (size(lambda_minus) /= size(n_minus_updated))) then
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_neutral) /= size(lambda_minus)"
       print *, "or size(lambda_minus) == size(n_minus_updated)"
       print *, "...aborting...", size(lambda_neutral), size(lambda_minus), size(n_minus_updated)
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.

    ! F = integrate_phi_spherical(exp(lambda_neutral))

    ! CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    ! CF3OO = integrate_phi_spherical(exp(lambda_neutral) * CF3 * F * F) !Note: F = "O" as both the F and O atoms and modelled as neutral

    ! n_minus_updated = bulk_density * exp(lambda_minus) * CF3OO * CF3OO

    F = integrate_phi_spherical(exp(lambda_neutral))
    O = integrate_phi_spherical(exp(lambda_minus))

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F * F * F)

    COOCF3 = integrate_phi_spherical(exp(lambda_neutral) * O * O * CF3)

    NCOOCF3 = integrate_phi_spherical(exp(lambda_neutral) * COOCF3)

    NCOOCF3COCF3 = integrate_phi_spherical(exp(lambda_neutral) * CF3 * O * NCOOCF3)

    n_minus_updated = 4.0_dp * bulk_density * exp(lambda_minus) * NCOOCF3COCF3

    !Explicitly ensure symmetry
    do ij = 1, (size(n_minus_updated) - 1)/2
       n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    end do

    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdateC4MIMTFSINegativeBeadDensities_model2

  !****************************************************************
  !*****END OF UPDATE DENSITY ROUTINES
  !****************************************************************

  !****************************************************************
  !*****START OF IDEAL CHAIN TERM ROUTINES
  !****************************************************************

  function calculate_single_neutral_sphere_ideal_chain_term(n_neutral)
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp) :: calculate_single_neutral_sphere_ideal_chain_term

    real(dp), dimension(size(n_neutral)) :: integrand

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_neutral))

    integrand(start_z_index:end_z_index) = n_neutral(start_z_index:end_z_index) * (log(n_neutral(start_z_index:end_z_index)) - 1.0_dp)

    calculate_single_neutral_sphere_ideal_chain_term = integrate_z_cylindrical(integrand, unity_function) / beta

  end function calculate_single_neutral_sphere_ideal_chain_term


  function calculate_neutral_dimers_ideal_chain_term(lambda_neutral)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp) :: calculate_neutral_dimers_ideal_chain_term

    real(dp), dimension(size(lambda_neutral)) :: integrand
    real(dp), dimension(size(lambda_neutral)) :: integrand_with_lambda

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    integrand_with_lambda(:) = 0.0_dp

    ! integrand(:) = 4.0_dp * pi * (hs_diameter**2) * integrate_phi_spherical(exp(lambda_neutral))
    ! integrand_with_lambda(:) = 4.0_dp * pi * (hs_diameter**2) * integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    integrand(:) = integrate_phi_spherical(exp(lambda_neutral))
    integrand_with_lambda(:) = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    calculate_neutral_dimers_ideal_chain_term = bulk_density * integrate_z_cylindrical(&
         (integrand_with_lambda(:) + (integrand(:) * (lambda_neutral + log(bulk_density) - 1.0_dp))) * exp(lambda_neutral), unity_function ) / beta

  end function calculate_neutral_dimers_ideal_chain_term


  function calculate_PositiveMinusSpheres_ideal_chain_term(n_plus, n_minus)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_minus

    real(dp) :: calculate_PositiveMinusSpheres_ideal_chain_term

    real(dp), dimension(size(n_plus)) :: integrand_plus, integrand_minus

    integer :: start_z_index
    integer :: end_z_index

    integrand_plus(:) = 0.0_dp
    integrand_minus(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_plus))

    integrand_plus(start_z_index:end_z_index) = n_plus(start_z_index:end_z_index) * (log(n_plus(start_z_index:end_z_index)) - 1.0_dp)
    integrand_minus(start_z_index:end_z_index) = n_minus(start_z_index:end_z_index) * (log(n_minus(start_z_index:end_z_index)) - 1.0_dp)

    calculate_PositiveMinusSpheres_ideal_chain_term = integrate_z_cylindrical(integrand_plus, unity_function) / beta +&
         integrate_z_cylindrical(integrand_minus, unity_function) / beta

  end function calculate_PositiveMinusSpheres_ideal_chain_term

  function calculate_PositiveNeutralMinusSpheres_ideal_chain_term(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus

    real(dp) :: calculate_PositiveNeutralMinusSpheres_ideal_chain_term

    real(dp), dimension(size(n_plus)) :: integrand_plus, integrand_neutral, integrand_minus

    integer :: start_z_index
    integer :: end_z_index

    integrand_plus(:) = 0.0_dp
    integrand_neutral(:) = 0.0_dp
    integrand_minus(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_plus))

    integrand_plus(start_z_index:end_z_index) = n_plus(start_z_index:end_z_index) * (log(n_plus(start_z_index:end_z_index)) - 1.0_dp)
    integrand_neutral(start_z_index:end_z_index) = n_neutral(start_z_index:end_z_index) * (log(n_neutral(start_z_index:end_z_index)) - 1.0_dp)
    integrand_minus(start_z_index:end_z_index) = n_minus(start_z_index:end_z_index) * (log(n_minus(start_z_index:end_z_index)) - 1.0_dp)

    calculate_PositiveNeutralMinusSpheres_ideal_chain_term = integrate_z_cylindrical(integrand_plus, unity_function) / beta +&
         integrate_z_cylindrical(integrand_neutral, unity_function) / beta + integrate_z_cylindrical(integrand_minus, unity_function) / beta

  end function calculate_PositiveNeutralMinusSpheres_ideal_chain_term

  function calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus, n_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(in) :: n_minus
    real(dp), intent(in) :: Donnan_potential
    real(dp) :: calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term

    real(dp), dimension(size(lambda_neutral)) :: integrand
    real(dp), dimension(size(lambda_neutral)) :: integrand_with_lambda

    integer :: start_z_index
    integer :: end_z_index

    real(dp) :: anion_contribution

    integrand(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_minus))

    !print *, "n_minus = ", n_minus
    !print *, "lambda = ", bulk_density*exp(lambda_minus)*exp(beta*negative_oligomer_charge*Donnan_potential)
    !call abort()
    !n_minus(start_z_index:end_z_index) * (log(n_minus(start_z_index:end_z_index)) - 1.0_dp)

    integrand(:) = (bulk_density*exp(lambda_minus)*exp(beta*negative_oligomer_charge*Donnan_potential)) &
         *(log(bulk_density*exp(lambda_minus)*exp(beta*negative_oligomer_charge*Donnan_potential)) - 1.0_dp)/beta

    anion_contribution = integrate_z_cylindrical(integrand, unity_function)

    integrand(:) = integrate_phi_spherical(exp(lambda_plus))
    integrand_with_lambda(:) = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term = anion_contribution &
         + (bulk_density/beta * exp(beta*positive_oligomer_charge*Donnan_potential) * integrate_z_cylindrical(integrand_with_lambda(:)*exp(lambda_neutral), unity_function)) &
         + (bulk_density/beta * exp(beta*positive_oligomer_charge*Donnan_potential) * integrate_z_cylindrical(integrand(:)*lambda_neutral(:)*exp(lambda_neutral), unity_function)) &
         + (bulk_density/beta * exp(beta*positive_oligomer_charge*Donnan_potential) * integrate_z_cylindrical(integrand(:)*log(bulk_density)*exp(lambda_neutral), unity_function)) &
         + (bulk_density/beta * exp(beta*positive_oligomer_charge*Donnan_potential) * integrate_z_cylindrical(integrand(:)*beta*positive_oligomer_charge*Donnan_potential*exp(lambda_neutral), unity_function)) &
         - (bulk_density/beta * exp(beta*positive_oligomer_charge*Donnan_potential) * integrate_z_cylindrical(integrand(:)*exp(lambda_neutral), unity_function))

    ! (bulk_density * exp(beta*positive_oligomer_charge*Donnan_potential) * &
    !      integrate_z_cylindrical((integrand_with_lambda(:) + (integrand(:) * (lambda_neutral + log(bulk_density) - 1.0_dp + beta*positive_oligomer_charge*Donnan_potential))) * exp(lambda_neutral), unity_function) / beta)

    ! calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term = anion_contribution + (bulk_density * &
    !      integrate_z_cylindrical((integrand_with_lambda(:) + (integrand(:) * (lambda_neutral + log(bulk_density) - 1.0_dp))) * exp(lambda_neutral), unity_function) / beta)

  end function calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term


  function calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp) :: calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term

    real(dp), dimension(size(lambda_neutral)) :: c1, c2, c3, c4, a1, a2
    real(dp), dimension(size(lambda_neutral)) :: c1_lambda, c2_lambda, c3_lambda, c4_lambda, a1_lambda, a2_lambda
    real(dp), dimension(size(lambda_neutral)) :: cation_integrand, anion_integrand

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    cation_integrand  = 0.0_dp
    anion_integrand = 0.0_dp

    c4 = integrate_phi_spherical(exp(lambda_neutral))
    c4_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c3 = integrate_phi_spherical(exp(lambda_neutral) * c4)
    c3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c4_lambda + c4*lambda_neutral))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_plus
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_plus) * cation_integrand, unity_function)*exp(beta*Donnan_potential*positive_oligomer_charge)


    a1 = integrate_phi_spherical(exp(lambda_minus))
    a1_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = a1_lambda + a1*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + a1*lambda_minus
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)*exp(beta*Donnan_potential*negative_oligomer_charge)

    calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term = anion_contribution + cation_contribution

  end function calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term


  function calculate_C4MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C4MIMBF4_ideal_chain_term

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    integer :: array_size

    real(dp), dimension(:), allocatable :: c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(:), allocatable :: c910, c4p
    real(dp), dimension(:), allocatable :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(:), allocatable :: c910_lambda, c4p_lambda

    real(dp), dimension(:), allocatable :: a1234
    real(dp), dimension(:), allocatable :: a1234_lambda

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) == size(lambda_neutral)) .or. (size(lambda_plus) == size(lambda_minus))) then
       allocate(c8(array_size), c7(array_size), c6(array_size), c5(array_size), c4(array_size), c3(array_size), c2(array_size), &
            c910(array_size), c4p(array_size), c8_lambda(array_size), c7_lambda(array_size), c6_lambda(array_size), c5_lambda(array_size), &
            c4_lambda(array_size), c3_lambda(array_size), c2_lambda(array_size), c910_lambda(array_size), c4p_lambda(array_size), &
            a1234(array_size), a1234_lambda(array_size), anion_integrand(array_size), cation_integrand(array_size))
    else
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    a1234 = integrate_phi_spherical(exp(lambda_minus))
    a1234_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = ( 4.0_dp*(a1234**3) * a1234_lambda ) + (a1234**4)*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + ((a1234**4)*lambda_minus)
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)*exp(beta*Donnan_potential*negative_oligomer_charge)

    c8 = integrate_phi_spherical(exp(lambda_neutral))
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C4MIMBF4_ideal_chain_term =  cation_contribution + anion_contribution
    !print *, "calculate_C4MIMBF4_ideal_chain_term = ", calculate_C4MIMBF4_ideal_chain_term
    !call abort()

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c8, c7, c6, c5, c4, c3, c2, c910, c4p, c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, &
         c3_lambda, c2_lambda, c910_lambda, c4p_lambda, a1234, a1234_lambda, anion_integrand, cation_integrand)

  end function calculate_C4MIMBF4_ideal_chain_term

  function calculate_C2MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C2MIMBF4_ideal_chain_term

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda

    real(dp), dimension(:), allocatable :: a1234
    real(dp), dimension(:), allocatable :: a1234_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if


    a1234 = integrate_phi_spherical(exp(lambda_minus))
    a1234_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = ( 4.0_dp*(a1234**3) * a1234_lambda ) + (a1234**4)*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + ((a1234**4)*lambda_minus)
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)*exp(beta*Donnan_potential*negative_oligomer_charge)


    c6 = integrate_phi_spherical(exp(lambda_neutral))
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C2MIMBF4_ideal_chain_term = anion_contribution + cation_contribution

  end function calculate_C2MIMBF4_ideal_chain_term

  
  
  function calculate_C6MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C6MIMBF4_ideal_chain_term

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda

    real(dp), dimension(:), allocatable :: a1234
    real(dp), dimension(:), allocatable :: a1234_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if


    a1234 = integrate_phi_spherical(exp(lambda_minus))
    a1234_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = ( 4.0_dp*(a1234**3) * a1234_lambda ) + (a1234**4)*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + ((a1234**4)*lambda_minus)
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)*exp(beta*Donnan_potential*negative_oligomer_charge)

    c12 = integrate_phi_spherical(exp(lambda_neutral))
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C6MIMBF4_ideal_chain_term = anion_contribution + cation_contribution

  end function calculate_C6MIMBF4_ideal_chain_term


  function calculate_C8MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C8MIMBF4_ideal_chain_term

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c13, c14, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda, c13_lambda, c14_lambda

    real(dp), dimension(:), allocatable :: a1234
    real(dp), dimension(:), allocatable :: a1234_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    a1234 = integrate_phi_spherical(exp(lambda_minus))
    a1234_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = ( 4.0_dp*(a1234**3) * a1234_lambda ) + (a1234**4)*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + ((a1234**4)*lambda_minus)
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)*exp(beta*Donnan_potential*negative_oligomer_charge)


    c14 = integrate_phi_spherical(exp(lambda_neutral))
    c14_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)
    c13_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c14_lambda + c14*lambda_neutral))

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c13_lambda + c13*lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C8MIMBF4_ideal_chain_term = anion_contribution + cation_contribution

  end function calculate_C8MIMBF4_ideal_chain_term


  function calculate_C10MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C10MIMBF4_ideal_chain_term

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c13, c14, c15, c16, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda, c13_lambda, c14_lambda, c15_lambda, c16_lambda

    real(dp), dimension(:), allocatable :: a1234
    real(dp), dimension(:), allocatable :: a1234_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if


    a1234 = integrate_phi_spherical(exp(lambda_minus))
    a1234_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = ( 4.0_dp*(a1234**3) * a1234_lambda ) + (a1234**4)*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + ((a1234**4)*lambda_minus)
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)*exp(beta*Donnan_potential*negative_oligomer_charge)

    c16 = integrate_phi_spherical(exp(lambda_neutral))
    c16_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16)
    c15_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c16_lambda + c16*lambda_neutral))

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)
    c14_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c15_lambda + c15*lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)
    c13_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c14_lambda + c14*lambda_neutral))

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c13_lambda + c13*lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C10MIMBF4_ideal_chain_term = anion_contribution + cation_contribution

  end function calculate_C10MIMBF4_ideal_chain_term















  
  function calculate_C4MIMTFSI_ideal_chain_term_model1(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C4MIMTFSI_ideal_chain_term_model1

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = F
    O_lambda = F_lambda

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_minus) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_minus) * (CF3COO_lambda + CF3COO*lambda_minus))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)


    c8 = integrate_phi_spherical(exp(lambda_neutral))
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C4MIMTFSI_ideal_chain_term_model1 = anion_contribution + cation_contribution

  end function calculate_C4MIMTFSI_ideal_chain_term_model1


  function calculate_C6MIMTFSI_ideal_chain_term_model1(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C6MIMTFSI_ideal_chain_term_model1

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = F
    O_lambda = F_lambda

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_minus) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_minus) * (CF3COO_lambda + CF3COO*lambda_minus))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c12 = integrate_phi_spherical(exp(lambda_neutral))
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C6MIMTFSI_ideal_chain_term_model1 = anion_contribution + cation_contribution

  end function calculate_C6MIMTFSI_ideal_chain_term_model1


  function calculate_C8MIMTFSI_ideal_chain_term_model1(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C8MIMTFSI_ideal_chain_term_model1

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c13, c14, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda, c13_lambda, c14_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = F
    O_lambda = F_lambda

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_minus) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_minus) * (CF3COO_lambda + CF3COO*lambda_minus))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c14 = integrate_phi_spherical(exp(lambda_neutral))
    c14_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)
    c13_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c14_lambda + c14*lambda_neutral))

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c13_lambda + c13*lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C8MIMTFSI_ideal_chain_term_model1 = anion_contribution + cation_contribution

  end function calculate_C8MIMTFSI_ideal_chain_term_model1


  function calculate_C10MIMTFSI_ideal_chain_term_model1(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C10MIMTFSI_ideal_chain_term_model1

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c13, c14, c15, c16, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda, c13_lambda, c14_lambda, c15_lambda, c16_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = F
    O_lambda = F_lambda

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_minus) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_minus) * (CF3COO_lambda + CF3COO*lambda_minus))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c16 = integrate_phi_spherical(exp(lambda_neutral))
    c16_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16)
    c15_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c16_lambda + c16*lambda_neutral))

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)
    c14_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c15_lambda + c15*lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)
    c13_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c14_lambda + c14*lambda_neutral))

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c13_lambda + c13*lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C10MIMTFSI_ideal_chain_term_model1 = anion_contribution + cation_contribution

  end function calculate_C10MIMTFSI_ideal_chain_term_model1


  function calculate_C2MIMTFSI_ideal_chain_term_model1(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C2MIMTFSI_ideal_chain_term_model1

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = F
    O_lambda = F_lambda

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_minus) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_minus) * (CF3COO_lambda + CF3COO*lambda_minus))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)


    c6 = integrate_phi_spherical(exp(lambda_neutral))
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C2MIMTFSI_ideal_chain_term_model1 = anion_contribution + cation_contribution

  end function calculate_C2MIMTFSI_ideal_chain_term_model1


  function calculate_C4MIMTFSI_ideal_chain_term_model2(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C4MIMTFSI_ideal_chain_term_model2

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda


    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = integrate_phi_spherical(exp(lambda_minus))
    O_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_neutral) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COO_lambda + CF3COO*lambda_neutral))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)


    c8 = integrate_phi_spherical(exp(lambda_neutral))
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C4MIMTFSI_ideal_chain_term_model2 = anion_contribution + cation_contribution

  end function calculate_C4MIMTFSI_ideal_chain_term_model2

  
  function calculate_C6MIMTFSI_ideal_chain_term_model2(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C6MIMTFSI_ideal_chain_term_model2

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = integrate_phi_spherical(exp(lambda_minus))
    O_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_neutral) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COO_lambda + CF3COO*lambda_neutral))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c12 = integrate_phi_spherical(exp(lambda_neutral))
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C6MIMTFSI_ideal_chain_term_model2 = anion_contribution + cation_contribution

  end function calculate_C6MIMTFSI_ideal_chain_term_model2


  function calculate_C8MIMTFSI_ideal_chain_term_model2(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C8MIMTFSI_ideal_chain_term_model2

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c13, c14, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda, c13_lambda, c14_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = integrate_phi_spherical(exp(lambda_minus))
    O_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_neutral) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COO_lambda + CF3COO*lambda_neutral))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c14 = integrate_phi_spherical(exp(lambda_neutral))
    c14_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)
    c13_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c14_lambda + c14*lambda_neutral))

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c13_lambda + c13*lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C8MIMTFSI_ideal_chain_term_model2 = anion_contribution + cation_contribution

  end function calculate_C8MIMTFSI_ideal_chain_term_model2


  function calculate_C10MIMTFSI_ideal_chain_term_model2(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C10MIMTFSI_ideal_chain_term_model2

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c11, c12, c13, c14, c15, c16, c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda, c11_lambda, c12_lambda, c13_lambda, c14_lambda, c15_lambda, c16_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = integrate_phi_spherical(exp(lambda_minus))
    O_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_neutral) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COO_lambda + CF3COO*lambda_neutral))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c16 = integrate_phi_spherical(exp(lambda_neutral))
    c16_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c15 = integrate_phi_spherical(exp(lambda_neutral) * c16)
    c15_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c16_lambda + c16*lambda_neutral))

    c14 = integrate_phi_spherical(exp(lambda_neutral) * c15)
    c14_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c15_lambda + c15*lambda_neutral))

    c13 = integrate_phi_spherical(exp(lambda_neutral) * c14)
    c13_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c14_lambda + c14*lambda_neutral))

    c12 = integrate_phi_spherical(exp(lambda_neutral) * c13)
    c12_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c13_lambda + c13*lambda_neutral))

    c11 = integrate_phi_spherical(exp(lambda_neutral) * c12)
    c11_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c12_lambda + c12*lambda_neutral))

    c8 = integrate_phi_spherical(exp(lambda_neutral) * c11)
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c11_lambda + c11*lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C10MIMTFSI_ideal_chain_term_model2 = anion_contribution + cation_contribution

  end function calculate_C10MIMTFSI_ideal_chain_term_model2


  function calculate_C2MIMTFSI_ideal_chain_term_model2(lambda_plus, lambda_neutral, lambda_minus, Donnan_potential)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_C2MIMTFSI_ideal_chain_term_model2

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    integer :: array_size

    real(dp), dimension(size(lambda_plus)) :: c6, c5, c4, c3, c2
    real(dp), dimension(size(lambda_plus)) :: c910, c4p
    real(dp), dimension(size(lambda_plus)) :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(size(lambda_plus)) :: c910_lambda, c4p_lambda

    real(dp), dimension(size(lambda_plus)) :: F, F3, O, CF3, CF3OO, CF3COO, two_CF3COO
    real(dp), dimension(size(lambda_plus)) :: F_lambda, F3_lambda, O_lambda, CF3_lambda, CF3OO_lambda, CF3COO_lambda, two_CF3COO_lambda

    real(dp), dimension(size(lambda_plus)) :: CF3COON, CF3COON_lambda, CF3COONOO, CF3COONOO_lambda, CF3COONCOO, CF3COONCOO_lambda
    real(dp), dimension(size(lambda_plus)) ::CF3COONCOOF2, CF3COONCOOF2_lambda, CF3COONCOOCF2, CF3COONCOOCF2_lambda

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    F = integrate_phi_spherical(exp(lambda_neutral))
    F_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    O = integrate_phi_spherical(exp(lambda_minus))
    O_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    F3 = F ** 3
    F3_lambda = 3.0_dp * (F_lambda * F * F)

    CF3 = integrate_phi_spherical(exp(lambda_neutral) * F3)
    CF3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (F3_lambda + F3*lambda_neutral))

    CF3OO = O * O * CF3
    CF3OO_lambda = (2.0_dp * (O_lambda * O * CF3)) + (CF3_lambda * O * O)

    CF3COO = integrate_phi_spherical(exp(lambda_neutral) * CF3OO)
    CF3COO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3OO_lambda + CF3OO*lambda_neutral))

    ! two_CF3COO = CF3COO * CF3COO
    ! two_CF3COO_lambda = 2.0_dp * (CF3COO * CF3COO_lambda)

    ! anion_integrand = two_CF3COO_lambda + (two_CF3COO * (log(bulk_density) - 1.0_dp)) + (two_CF3COO*exp(lambda_minus))
    ! anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    CF3COON = integrate_phi_spherical(exp(lambda_neutral) * CF3COO)
    CF3COON_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COO_lambda + CF3COO*lambda_neutral))

    CF3COONOO = CF3COON * O * O
    CF3COONOO_lambda = (CF3COON_lambda * O * O) + (2.0_dp * O_lambda * O * CF3COON)

    CF3COONCOO = integrate_phi_spherical(exp(lambda_neutral) * CF3COONOO)
    CF3COONCOO_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONOO_lambda + CF3COONOO*lambda_neutral))

    CF3COONCOOF2 = F * F * CF3COONCOO
    CF3COONCOOF2_lambda = (2.0_dp * F_lambda * F * CF3COONCOO) + (CF3COONCOO_lambda * F * F)

    CF3COONCOOCF2 = integrate_phi_spherical(exp(lambda_neutral) * CF3COONCOOF2)
    CF3COONCOOCF2_lambda = integrate_phi_spherical(exp(lambda_neutral) * (CF3COONCOOF2_lambda + CF3COONCOOF2*lambda_neutral))

    anion_integrand = CF3COONCOOCF2_lambda + CF3COONCOOCF2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*negative_oligomer_charge)) + CF3COONCOOCF2*lambda_neutral
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * anion_integrand, unity_function) * exp(beta*Donnan_potential*negative_oligomer_charge)

    c6 = integrate_phi_spherical(exp(lambda_neutral))
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp + (beta*Donnan_potential*positive_oligomer_charge)) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function) * exp(beta*Donnan_potential*positive_oligomer_charge)

    calculate_C2MIMTFSI_ideal_chain_term_model2 = anion_contribution + cation_contribution

  end function calculate_C2MIMTFSI_ideal_chain_term_model2

  !****************************************************************
  !*****END OF IDEAL CHAIN TERM ROUTINES
  !****************************************************************

  !****************************************************************
  !*****START OF CHEMICAL POTENTIAL ROUTINES
  !****************************************************************

  function calculate_chem_potential_term_neutral_spheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_term_neutral_spheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus
    
    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend
    
    real(dp) :: lambda

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_term_neutral_spheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chem_potential_term_neutral_spheres = (1.0_dp/beta) * (log(bulk_density_neutral_beads) + lambda) * &
         integrate_z_cylindrical(n_neutral, unity_function)

  end function calculate_chem_potential_term_neutral_spheres


  function calculate_chem_potential_term_neutral_dimers(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    !real(dp), dimension(size(n_plus)) :: lambda_hs_end
    !real(dp), dimension(size(n_plus)) :: lambda_hs_nonend

    real(dp) :: calculate_chem_potential_term_neutral_dimers

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_term_neutral_spheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    !call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation, lambda_hs_end, lambda_hs_nonend)

    ! calculate_chem_potential_term_neutral_dimers = (1.0_dp/beta) * (log(bulk_density) + (2.0_dp * lambda)) * bulk_density * &
    !      integrate_z_cylindrical(integrate_phi_spherical(exp(lambda_neutral)) * exp(lambda_neutral), unity_function)

    calculate_chem_potential_term_neutral_dimers = (1.0_dp/beta) * (log(bulk_density) + (2.0_dp * lambda)) * &
         integrate_z_cylindrical(0.5_dp * n_neutral, unity_function)

  end function calculate_chem_potential_term_neutral_dimers


  function calculate_chem_potential_PositiveMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_PositiveMinusSpheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda_plus_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chem_potential_PositiveMinusSpheres = &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk) * &
         integrate_z_cylindrical(n_plus, unity_function)) + &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk) * &
         integrate_z_cylindrical(n_minus, unity_function))

  end function calculate_chem_potential_PositiveMinusSpheres

  function calculate_chem_potential_PositiveNeutralMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation, Donnan_potential)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_chem_potential_PositiveNeutralMinusSpheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chem_potential_PositiveNeutralMinusSpheres = &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + (Donnan_potential*positive_oligomer_charge)) * &
         integrate_z_cylindrical(n_plus, unity_function)) + &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_neutral_bulk) * &
         integrate_z_cylindrical(n_neutral, unity_function)) + &         
         ((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk + (Donnan_potential*negative_oligomer_charge)) * &
         integrate_z_cylindrical(n_minus, unity_function))

  end function calculate_chem_potential_PositiveNeutralMinusSpheres


  function calculate_chem_potential_PositiveNeutralDimerMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation, Donnan_potential)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation
    real(dp), intent(in) :: Donnan_potential

    !real(dp), dimension(size(n_plus)) :: lambda_hs_end
    !real(dp), dimension(size(n_plus)) :: lambda_hs_nonend

    real(dp) :: calculate_chem_potential_PositiveNeutralDimerMinusSpheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp), dimension(size(n_neutral)) :: integrand, integrand_with_lambda, c1, c2

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk
    real(dp) :: positive_oligomer_density

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    ! calculate_chem_potential_PositiveNeutralDimerMinusSpheres = &
    !      (((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + lambda_neutral_bulk) + (Donnan_potential*positive_oligomer_charge)) * &
    !      integrate_z_cylindrical((n_plus + n_neutral)/2.0_dp, unity_function)) + &
    !      (((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk) + (Donnan_potential*negative_oligomer_charge)) * &
    !      integrate_z_cylindrical(n_minus, unity_function))


    !call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation, lambda_hs_end, lambda_hs_nonend)


    !integrand(:) = integrate_phi_spherical(exp(lambda_neutral))
    !integrand_with_lambda(:) = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    !positive_oligomer_density = (bulk_density * integrate_z_cylindrical(&
    !     (integrand_with_lambda(:) + (integrand(:) * (lambda_plus + log(bulk_density) - 1.0_dp))) * exp(lambda_plus), unity_function ) / beta)
    !*exp(beta*positive_oligomer_density*Donnan_potential)

    ! calculate_chem_potential_PositiveNeutralDimerMinusSpheres = &
    !      (((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + lambda_neutral_bulk) + (Donnan_potential*positive_oligomer_charge)) * &
    !      positive_oligomer_density) + &
    !      (((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk) + (Donnan_potential*negative_oligomer_charge)) * &
    !      integrate_z_cylindrical(n_minus, unity_function))

    !calculate_chem_potential_PositiveNeutralDimerMinusSpheres = &
    !     (((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + lambda_neutral_bulk)) * &
    !     positive_oligomer_density) + &
    !     (((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk)) * &
    !     integrate_z_cylindrical(n_minus, unity_function))

    c1 = integrate_phi_spherical(exp(lambda_plus))
    c2 = bulk_density * exp(lambda_neutral) * c1

    calculate_chem_potential_PositiveNeutralDimerMinusSpheres = &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + lambda_neutral_bulk + (Donnan_potential * positive_oligomer_charge)) *&
         integrate_z_cylindrical(c2*exp(beta*Donnan_potential*positive_oligomer_charge), unity_function)) + &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk + (Donnan_potential * negative_oligomer_charge))*&
         integrate_z_cylindrical(bulk_density*exp(lambda_minus)*exp(beta*Donnan_potential*negative_oligomer_charge), unity_function))


    ! calculate_chem_potential_PositiveNeutralDimerMinusSpheres = &
    !      ((((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + lambda_neutral_bulk))) *&
    !      integrate_z_cylindrical(c2, unity_function) ) + &
    !      ((((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk)))*&
    !      integrate_z_cylindrical(bulk_density*exp(lambda_minus), unity_function))



  end function calculate_chem_potential_PositiveNeutralDimerMinusSpheres


  function calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    !call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation, lambda_hs_end, lambda_hs_nonend)

    ! calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * (&
    !      (log(bulk_density) + ((2.0_dp * lambda_plus_bulk) + (2.0_dp * lambda_neutral_bulk))) * &
    !      (integrate_z_cylindrical((n_plus + n_neutral)/4.0_dp, unity_function)) + &
    !      (integrate_z_cylindrical(n_minus/2.0_dp, unity_function) * (log(bulk_density) + 2.0_dp * lambda_minus_bulk)))

    calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * (&
         (log(bulk_density)*integrate_z_cylindrical(n_plus/2.0_dp, unity_function)) + &
         (log(bulk_density)*integrate_z_cylindrical(n_minus/2.0_dp, unity_function)) + &
         (lambda_plus_bulk*integrate_z_cylindrical(n_plus, unity_function)) + &
         (lambda_neutral_bulk*integrate_z_cylindrical(n_neutral, unity_function)) + &
         (lambda_minus_bulk*integrate_z_cylindrical(n_minus, unity_function)) + &
         ((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_plus/2.0_dp, unity_function)) + &
         ((negative_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_minus/2.0_dp, unity_function)) &
         )

    ! calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * (&
    !      (log(bulk_density) + ((2.0_dp * lambda_plus_bulk) + (2.0_dp * lambda_neutral_bulk))) * &
    !      (integrate_z_cylindrical((n_plus + n_neutral)/4.0_dp, unity_function)))

    ! calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * &
    !      (log(bulk_density) + ((2.0_dp * lambda_minus_bulk))) * &
    !      (integrate_z_cylindrical(0.5_dp * n_minus, unity_function))

  end function calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer


  function calculate_chem_potential_C4MIMBF4(n_plus, n_neutral, n_minus, ith_plate_separation, Donnan_potential)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation
    real(dp), intent(in) :: Donnan_potential

    real(dp) :: calculate_chem_potential_C4MIMBF4

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    !call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation, lambda_hs_end, lambda_hs_nonend)

    !calculate_chem_potential_C4MIMBF4 = 0.0_dp

    !calculate_chem_potential_C4MIMBF4 = (1.0_dp/beta) * (&
    !     (log(bulk_density) + ((5.0_dp * lambda_plus_bulk) + (5.0_dp * lambda_neutral_bulk) + (Donnan_potential*positive_oligomer_charge))) * &
    !     (integrate_z_cylindrical((n_plus + n_neutral)/10.0_dp, unity_function)) + &
    !     (integrate_z_cylindrical(n_minus/5.0_dp, unity_function) * (log(bulk_density) + (5.0_dp * lambda_minus_bulk) + (Donnan_potential*negative_oligomer_charge))))



    calculate_chem_potential_C4MIMBF4 = (1.0_dp/beta) * (&
         (log(bulk_density)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
         (log(bulk_density)*integrate_z_cylindrical(n_minus/5.0_dp, unity_function)) + &
         (lambda_plus_bulk*integrate_z_cylindrical(n_plus, unity_function)) + &
         (lambda_neutral_bulk*integrate_z_cylindrical(n_neutral, unity_function)) + &
         (lambda_minus_bulk*integrate_z_cylindrical(n_minus, unity_function)) + &
         ((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
         ((negative_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_minus/5.0_dp, unity_function)) &
         )




    ! calculate_chem_potential_C4MIMBF4 = ((log(bulk_density)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function))/beta) + &
    !      ((log(bulk_density)*integrate_z_cylindrical(n_minus, unity_function))/beta) + &
    !      (((lambda_plus_bulk)*integrate_z_cylindrical(n_plus, unity_function))/beta) + &
    !      ((lambda_neutral_bulk*integrate_z_cylindrical(n_neutral, unity_function))/beta) + &
    !      (((lambda_minus_bulk)*integrate_z_cylindrical(n_minus, unity_function))/beta)! + &

    !((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_plus, unity_function)) + &
    !     ((negative_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_minus, unity_function)) + &
    !     ((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_neutral, unity_function))

    ! print *, (1.0_dp/beta) * (&
    !      (log(bulk_density)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
    !      (log(bulk_density)*integrate_z_cylindrical(n_minus/5.0_dp, unity_function)) + &
    !      ((lambda_plus_bulk)*integrate_z_cylindrical(n_plus, unity_function)) + &
    !      (lambda_neutral_bulk*integrate_z_cylindrical(n_neutral, unity_function)) + &
    !      ((lambda_minus_bulk)*integrate_z_cylindrical(n_minus, unity_function)))


    ! print *, ((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function))
    ! print *, ((negative_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_minus/5.0_dp, unity_function))
    ! print *, calculate_chem_potential_C4MIMBF4

    ! print *, log(bulk_density), lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk
    ! print *, integrate_z_cylindrical(n_plus/5.0_dp, unity_function), integrate_z_cylindrical(n_minus/5.0_dp, unity_function)
    ! print *, Donnan_potential, positive_oligomer_charge, negative_oligomer_charge

    ! call abort()

  end function calculate_chem_potential_C4MIMBF4

  function calculate_chem_potential_C4MIMTFSI_model1(n_plus, n_neutral, n_minus, ith_plate_separation, Donnan_potential)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation
    real(dp), intent(in) :: Donnan_potential

    !real(dp), dimension(size(n_plus)) :: lambda_hs_end
    !real(dp), dimension(size(n_plus)) :: lambda_hs_nonend

    real(dp) :: calculate_chem_potential_C4MIMTFSI_model1

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    !call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation, lambda_hs_end, lambda_hs_nonend)

    !calculate_chem_potential_C4MIMBF4 = 0.0_dp

    ! calculate_chem_potential_C4MIMTFSI_model1 = (1.0_dp/beta) * (&
    !      ((log(bulk_density) + ((5.0_dp * lambda_plus_bulk) + (5.0_dp * lambda_neutral_bulk))) * &
    !      (integrate_z_cylindrical((n_plus + ((5.0_dp/19.0_dp)*n_neutral))/10.0_dp, unity_function))) + &
    !      (integrate_z_cylindrical((n_minus + ((14.0_dp/19.0_dp)*n_neutral))/15.0_dp, unity_function) * (log(bulk_density) + lambda_minus_bulk + (lambda_neutral_bulk*14.0_dp))))

    calculate_chem_potential_C4MIMTFSI_model1 = (1.0_dp/beta) * (&
         (log(bulk_density)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
         (log(bulk_density)*integrate_z_cylindrical(n_minus, unity_function)) + &
         (lambda_plus_bulk*integrate_z_cylindrical(n_plus, unity_function)) + &
         (lambda_neutral_bulk*integrate_z_cylindrical(n_neutral, unity_function)) + &
         (lambda_minus_bulk*integrate_z_cylindrical(n_minus, unity_function)) + &
         ((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
         ((negative_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_minus, unity_function)) &
         )


  end function calculate_chem_potential_C4MIMTFSI_model1


  function calculate_chem_potential_C4MIMTFSI_model2(n_plus, n_neutral, n_minus, ith_plate_separation, Donnan_potential)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation
    real(dp), intent(in) :: Donnan_potential

    !real(dp), dimension(size(n_plus)) :: lambda_hs_end
    !real(dp), dimension(size(n_plus)) :: lambda_hs_nonend

    real(dp) :: calculate_chem_potential_C4MIMTFSI_model2

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp), dimension(size(n_neutral)) :: lambda_hs_end
    real(dp), dimension(size(n_minus)) :: lambda_hs_nonend
    real(dp), dimension(size(n_neutral)) :: n_hs_end
    real(dp), dimension(size(n_minus)) :: n_hs_nonend

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end, n_hs_end, lambda_hs_nonend, n_hs_nonend, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    !call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation, lambda_hs_end, lambda_hs_nonend)

    !calculate_chem_potential_C4MIMBF4 = 0.0_dp

    ! calculate_chem_potential_C4MIMTFSI_model1 = (1.0_dp/beta) * (&
    !      ((log(bulk_density) + ((5.0_dp * lambda_plus_bulk) + (5.0_dp * lambda_neutral_bulk))) * &
    !      (integrate_z_cylindrical((n_plus + ((5.0_dp/19.0_dp)*n_neutral))/10.0_dp, unity_function))) + &
    !      (integrate_z_cylindrical((n_minus + ((14.0_dp/19.0_dp)*n_neutral))/15.0_dp, unity_function) * (log(bulk_density) + lambda_minus_bulk + (lambda_neutral_bulk*14.0_dp))))

    calculate_chem_potential_C4MIMTFSI_model2 = (1.0_dp/beta) * (&
         (log(bulk_density)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
         (log(bulk_density)*integrate_z_cylindrical(n_minus/4.0_dp, unity_function)) + &
         (lambda_plus_bulk*integrate_z_cylindrical(n_plus, unity_function)) + &
         (lambda_neutral_bulk*integrate_z_cylindrical(n_neutral, unity_function)) + &
         (lambda_minus_bulk*integrate_z_cylindrical(n_minus, unity_function)) + &
         ((positive_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_plus/5.0_dp, unity_function)) + &
         ((negative_oligomer_charge*Donnan_potential)*integrate_z_cylindrical(n_minus/4.0_dp, unity_function)) &
         )


  end function calculate_chem_potential_C4MIMTFSI_model2

  !****************************************************************
  !*****END OF CHEMICAL POTENTIAL ROUTINES
  !****************************************************************

end module constructoligomers
