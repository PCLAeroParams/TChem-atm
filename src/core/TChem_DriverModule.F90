module TChemDriver

!interface
!  subroutine initialize(arg_chemfile) bind(c, name="initialize")
!    use iso_c_binding
!    character(kind=c_char) :: arg_chemfile(*)
!  end subroutine initialize
!  subroutine finalize() bind(c, name="finalize")
!  end subroutine finalize
!  subroutine getAllStateVectorHost(state_vector) bind(c)
!    use iso_c_binding
!    type(c_ptr) :: state_vector
!  end subroutine
!  function TChem_getNumberOfSpecies() bind(c, name="TChem_getNumberOfSpecies")
!    use iso_c_binding
!    integer(kind=c_int) :: TChem_getNumberOfSpecies
!!  end function
!  subroutine TChem_getSpeciesName(index) bind(c, name="TChem_getSpeciesName")
!    use iso_c_binding
!    integer(kind=c_int), intent(in) :: index
!    integer(kind=c_int), intent(out) :: lw
!    character(C_CHAR), intent(out)  :: w(*)
!  end subroutine 
!  subroutine TChem_getStateVector(array) bind(c, name="TChem_getStateVector")
!    use iso_c_binding
!    real(kind=c_double) :: array(*)
!  end subroutine
!  subroutine TChem_setStateVector(array) bind(c, name="TChem_setStateVector")
!    use iso_c_binding
!    real(kind=c_double) :: array(*)
!  end subroutine
!end interface

contains

!  subroutine timestep
!    use iso_c_binding
!
!    print*, 'we are doing some chemistry'
!
!  end subroutine timestep
!
!  subroutine TChem_initialize(chemFile)
!
!    use iso_c_binding
!
!    character(kind=c_char,len=*), intent(in) :: chemFile
!
!    print*, 'in TChemDriver: ', chemFile
!
!    call initialize(chemFile//c_null_char)
!
!  end subroutine

end module TChemDriver
