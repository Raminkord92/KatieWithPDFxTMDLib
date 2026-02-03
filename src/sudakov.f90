module katie_sudakov
use katie_quadpack
use avh_trivinp
implicit none
  private
  public :: sudakov_init,sudakov_g,sudakov_q

  real(kind(1d0)),parameter :: NA = 3d0
  real(kind(1d0)),parameter :: Cf = (NA*NA-1)/(2*NA)
  real(kind(1d0)),parameter :: small=epsilon(1d0)
  real(kind(1d0)),parameter :: r1PI = 3.1415926535897932384626433832795d0
  real(kind(1d0)),parameter :: r2PI = 2*r1PI
  integer,save :: Nf = 5

  real(kind(1d0)),parameter :: kt_freeze=0.5d0
  real(kind(1d0)),parameter :: mu_min=1d0 ,mu_max=1d4
  real(kind(1d0)),save :: mu=-1d0
  type(trivinp_type),save :: grid_g ,grid_q
  
  character(256),parameter :: tableDir=&
![sudakovTableDir
!]sudakovTableDir

contains

  subroutine test_it
  real(kind(1d0)) :: mu_in,kt_in
  integer :: npnt,ii
  character(256) :: val
  !call get_command_argument( 1 ,val ); read(val,*) kt_in
  !call get_command_argument( 2 ,val ); read(val,*) mu_in

  call sudakov_init( 'CT10nlo' ,5 )
  end subroutine

  
  subroutine sudakov_init( lhaSet ,Nflavors )
  character(*),intent(in) :: lhaSet
  integer,intent(in) :: Nflavors
  character(256) :: prefix
  integer :: ios
  Nf = Nflavors
  call InitPDFsetByName(trim(lhaSet)//'.LHgrid')
  prefix = trim(tableDir)//'sudakov_'//trim(lhaSet)
  open(33,file=trim(prefix)//'_g.dat',status='old',iostat=ios);close(33)
  if (ios.eq.0) then
    call grid_g%read_file(trim(prefix)//'_g.dat',33)
    call grid_q%read_file(trim(prefix)//'_q.dat',33)
  else
    write(*,'(A)') 'Creating grids for Sudakov factors...'
    call grid_g%create( integral_g ,0d0,1d0,200 ,log(mu_min),log(mu_max),200 )
    call grid_q%create( integral_q ,0d0,1d0,200 ,log(mu_min),log(mu_max),200 )
    call grid_g%write_file(trim(prefix)//'_g.dat',33)
    call grid_q%write_file(trim(prefix)//'_q.dat',33)
    write(*,'(A)') '...done.'
  endif
  end subroutine

  function alphas(kt) result(rslt)
  real(kind(1d0)),intent(in) :: kt
  real(kind(1d0)) :: rslt,alphasPDF
  if (kt.lt.kt_freeze) then
    rslt = alphasPDF(kt_freeze)
  else
    rslt = alphasPDF(kt)
  endif
  end function
  
  function intd_split_g(xx) result(rslt)
! xx is supposed to be kt/mu
  real(kind(1d0)),intent(in) :: xx
  real(kind(1d0)) :: rslt
  if (xx.le.small) then
    rslt = -12*log(small)-11
  else
    rslt = -12*log(xx) - (1-xx)*(11+xx*(20+xx*11))/(1+xx)**3
  endif
  rslt = (Cf/6)*rslt + Nf/3
  end function

  function intd_split_q(xx) result(rslt)
! xx is supposed to be kt/mu
  real(kind(1d0)),intent(in) :: xx
  real(kind(1d0)) :: rslt,hh
  if (xx.le.small) then
    rslt = -4*log(small)-3
  else
    rslt = 4*log(1+1/xx) - (3+2*xx)/(1+xx)**2 
  endif
  rslt = (Cf/2)*rslt
  end function

  function integrand_g(logx) result(rslt)
! logx is supposed to be log(kt/mu)
  real(kind(1d0)),intent(in) :: logx
  real(kind(1d0)) :: rslt,xx
  xx = exp(logx)
  rslt = alphas(mu*xx)*intd_split_g(xx)
  end function
  
  function integrand_q(logx) result(rslt)
! logx is supposed to be log(kt/mu)
  real(kind(1d0)),intent(in) :: logx
  real(kind(1d0)) :: rslt,xx
  xx = exp(logx)
  rslt = alphas(mu*xx)*intd_split_q(xx)
  end function

  function integral_g( kt_over_mu ,logmu ) result(rslt)
  real(kind(1d0)),intent(in) :: kt_over_mu,logmu
  real(kind(1d0)) :: rslt,logmin,logmax
  mu = exp(logmu)
  if (kt_over_mu.lt.small) then
    logmin = log(small)
  else
    logmin = log(kt_over_mu)
  endif
  logmax = 0
  rslt = integral( integrand_g ,logmin ,logmax )
  rslt = exp(-rslt/r1PI)
  end function
  
  function integral_q( kt_over_mu ,logmu ) result(rslt)
  real(kind(1d0)),intent(in) :: kt_over_mu,logmu
  real(kind(1d0)) :: rslt,logmin,logmax
  mu = exp(logmu)
  if (kt_over_mu.lt.small) then
    logmin = log(small)
  else
    logmin = log(kt_over_mu)
  endif
  logmax = 0
  rslt = integral( integrand_q ,logmin ,logmax )
  rslt = exp(-rslt/r1PI)
  end function
  
  function sudakov_g(kt,mu) result(rslt)
  real(kind(1d0)),intent(in) :: kt,mu
  real(kind(1d0)) :: rslt
  if (kt.gt.mu) then
    rslt = 1
  elseif (mu.lt.mu_min) then
    rslt = 0
  else
    rslt = grid_g%evaluate(kt/mu,log(mu))
  endif
  end function
  
  function sudakov_q(kt,mu) result(rslt)
  real(kind(1d0)),intent(in) :: kt,mu
  real(kind(1d0)) :: rslt
  if (kt.gt.mu) then
    rslt = 1
  elseif (mu.lt.mu_min) then
    rslt = 0
  else
    rslt = grid_q%evaluate(kt/mu,log(mu))
  endif
  end function

end module

!program test_sudakov
!use sudakov
!  call test_it
!end program
