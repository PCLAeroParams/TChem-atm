module TChemDriver

interface
  subroutine initialize_kokkos(arg_chemfile) bind(c, name="initialize_kokkos")
    use iso_c_binding
    character(kind=c_char) :: arg_chemfile(*)
  end subroutine initialize_kokkos
  subroutine finalize_kokkos() bind(c, name="finalize_kokkos")
  end subroutine finalize_kokkos
  subroutine getAllStateVectorHost(state_vector) bind(c)
    use iso_c_binding
    type(c_ptr) :: state_vector
  end subroutine
  function TChem_getNumberOfSpecies() bind(c, name="TChem_getNumberOfSpecies")
    use iso_c_binding
    integer(kind=c_int) :: TChem_getNumberOfSpecies
  end function
  function TChem_getSpeciesName(index) bind(c)
    use iso_c_binding
    integer(kind=c_int) :: index
    character(kind=c_char), dimension(50) :: TChem_TChem_getSpeciesName
  end function
  subroutine TChem_getStateVector(array) bind(c, name="TChem_getStateVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
  end subroutine
  subroutine TChem_setStateVector(array) bind(c, name="TChem_setStateVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
  end subroutine
     subroutine TChem_setAllStateVectorHost(arg_state_vector) bind(C, &
          name="TChem_setAllStateVectorHost")
       use iso_c_binding, only: c_ptr
       type(c_ptr), value :: arg_state_vector
     end subroutine TChem_setAllStateVectorHost
end interface

contains

  subroutine test
    use iso_c_binding

    print*, 'we are doing some chemistry'

  end subroutine test

  subroutine initialize(chemFile)

    use iso_c_binding

    character(kind=c_char,len=*), intent(in) :: chemFile

    print*, 'in TChemDriver: ', chemFile

    call initialize_kokkos(chemFile//c_null_char)

  end subroutine

end module TChemDriver
