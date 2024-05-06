module TChemDriver

interface
  subroutine normal_vec(n, x) bind(c, name="normal_vec")
    use iso_c_binding
    integer(kind=c_int), intent(in), value :: n
    real(kind=c_double)                    :: x(n)
  end subroutine normal_vec
  subroutine call_something(n) bind(c)
    use iso_c_binding
    integer(kind=c_int), intent(in) :: n
  end subroutine call_something

  subroutine do_something(n) bind(c)
    use iso_c_binding
    integer(kind=c_int), intent(in) :: n
  end subroutine do_something

  subroutine initialize_kokkos() bind(c)
!    use iso_c_binding
!    integer(kind=c_int) :: nSpec
  end subroutine initialize_kokkos
  subroutine finalize_kokkos() bind(c)
  end subroutine finalize_kokkos
  
end interface

contains

  subroutine test
    use iso_c_binding
    integer :: seed
    integer, parameter :: n = 5
    real(kind=c_double) :: x(n)

    print*, 'hello from TChem!'

    do seed=1,3
       print "(*(1x,f7.4))",x
    end do

    do seed=1,3
       call normal_vec(n, x)
       print "(*(1x,f7.4))",x
    end do

    call call_something(n)

    call do_something(n)

  end subroutine test

end module TChemDriver
