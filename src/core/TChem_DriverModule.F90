module TChemDriver

interface
  subroutine call_something(n) bind(c, name="call_something")
    use iso_c_binding
    integer(kind=c_int), intent(in) :: n
  end subroutine call_something
  subroutine do_something(n) bind(c, name="do_something")
    use iso_c_binding
    integer(kind=c_int), intent(in) :: n
  end subroutine do_something
  subroutine initialize_kokkos(arg_chemfile) bind(c, name="initialize_kokkos")
    use iso_c_binding
    character(kind=c_char) :: arg_chemfile(*)
  end subroutine initialize_kokkos
  subroutine finalize_kokkos() bind(c, name="finalize_kokkos")
  end subroutine finalize_kokkos
end interface

contains

  subroutine test
    use iso_c_binding
    integer, parameter :: n = 5

!    call call_something(n)
    call do_something(n)

  end subroutine test

  subroutine initialize(chemFile)

    use iso_c_binding

    character(kind=c_char,len=*), intent(in) :: chemFile

    print*, 'in TChemDriver: ', chemFile

    call initialize_kokkos(chemFile//c_null_char)

  end subroutine

  integer function TChem_getNumberOfSpecies()

     TChem_getNumberOfSpecies = 67

  end function

  subroutine TChem_getNamesOfSpecies()

  end subroutine

  subroutine TChem_getAllStateVectorHost() !state_vector)

    use iso_c_binding

    real(kind=c_double), dimension(:), allocatable, target :: state_vector

    integer,        parameter :: nbatch = 1
    integer :: len_state_vector

  end subroutine

end module TChemDriver
