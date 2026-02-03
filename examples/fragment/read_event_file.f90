!!=======================================================================
!!== List of protected variables, functions, and subroutines: ===========
!!=======================================================================
! integer :: Nsize            ! size parameter, max number of final-states
! integer :: Nfnst            ! total number of final states
! integer :: Ngroups          ! number of groups (>1 for multi-parton scattering)
! integer :: NfinalState(i)   ! number of final states in group i
! integer :: process(i)       ! id of the process of group i
! integer :: procNr           ! identical to process(1), convenient for SPS
! integer :: NfnstProc(process(i))  ! number of final-state particles in the
!                                   ! process of group i
! integer :: lhefInst(j,process(i)) ! particle id in LHEF encoding of initial-
!                                   ! state particle j the process of group i
! integer :: lhefFnst(j,process(i)) ! particle id in LHEF encoding of final-
!                                   ! state particle j the process of group i
! real(kind(1d0)) :: Ecm              ! total center-off-mass energy
! real(kind(1d0)) :: eventWeight      ! event weight
! real(kind(1d0)) :: xsection         ! total cross section
! real(kind(1d0)) :: errest           ! statistical error estimate on total cross section
! real(kind(1d0)) :: pInst(0:3,j,i)   ! initial-state momentum j of group i
! real(kind(1d0)) :: sInst(j,i)       ! initial-state squared momentum j of group i
! real(kind(1d0)) :: pFnst(0:3,j)     ! final-state momentum j
! real(kind(1d0)) :: sFnst(j)         ! final-state squared momentum j
! real(kind(1d0)) :: matrixElement(i) ! value of matrix element for group i
! real(kind(1d0)) :: partonLumin(i)   ! value of parton luminosity for group i
! real(kind(1d0)) :: alphaStrong(i)   ! value of strong coupling for group i
! real(kind(1d0)) :: muScale(i)       ! value energy scale for group i
! character(256)  :: filename         ! name (including path) of event file
! real(kind(1d0) function mass(p(0:3))            ! mass = sqrt(abs(square))
! real(kind(1d0) function pTrans(p(0:3))          ! transverse momentum
! real(kind(1d0) function ETrans(p(0:3))          ! transverse mass
! real(kind(1d0) function rapidity(p(0:3))        ! rapidity
! real(kind(1d0) function pseudoRap(p(0:3))       ! pseudo rapidity
! real(kind(1d0) function phi(p(0:3))             ! azimuthal angle
! real(kind(1d0) function deltaPhi(p(0:3),q(0:3)) ! difference between azim.angles
! real(kind(1d0) function deltaR(p(0:3),q(0:3))   ! delta R
! integer        :: unity(i)   ! the number i 
! character(1)   :: numeral(i) ! the number i as a character
! real(kind(1d0) :: r1PI       ! the number pi
! real(kind(1d0) :: r2PI       ! the number 2*pi
! subroutine sort_big2small(l,a,n) ! sort a(1:n) from big to small or small to big
! subroutine sort_small2big(l,a,n) ! the array l contains the associated permutation
!
! type :: histo_1d_type ! 1-dimensional histogram
! type :: histo_2d_type ! 2-dimensional histogram
!
! Initialize 1-dim histogram with:
!   call hst%init( left=0d0 ,right=2d0 ,Nbins=4 )
! or alternatively:
!   call hst%init( [0.0,0.5 ,0.5,1.0 ,1.0,1.5 ,1.5,2.0] )
! Data collection:
!   call hst%collect( value ,weight )
! Write histogram file:
!   call hst%write('somename')
!
! Initialize 2-dim histogram with:
!   call hst%init( xLow=0d0 ,xUpp=2d0 ,NxBins=4 ,yLow=-1d0 ,yUpp=1d0 ,Nybins=2 )
! or alternatively:
!   call hst%init( [0.0,0.5 ,0.5,1.0 ,1.0,1.5 ,1.5,2.0] ,[-1.0,0.0 ,0.0,1.0] )

!=======================================================================
program read_event_file ! ==================== DO NOT TOUCH THIS LINE ==
use read_events_mod ;use avh_random ! ======== DO NOT TOUCH THIS LINE ==
!=======================================================================
!== DECLARATIONS =======================================================
!=======================================================================

implicit none
type(histo_1d_type) :: h_pT(Nsize),h_eta(Nsize)
real(kind(1d0)) :: pT(Nsize),eta(Nsize),totalWeight
real(kind(1d0)) :: sumW0,sumW1,sumW2,xsection,errest
integer :: ii,pTordrd(Nsize),rapOrdrd(Nsize)

!=======================================================================
call open_file ! ============================= DO NOT TOUCH THIS LINE ==
!=======================================================================
!== INITIALIZATIONS ====================================================
!=======================================================================

call initmult() ! initialization of the subroutine swqmult
call rangen_init(seed=1234) ! initialization of random number generator

call h_pT(1)%init( left=0d0 ,right=200d0 ,Nbins=100 )
call h_pT(2)%init( left=0d0 ,right=200d0 ,Nbins=100 )
call h_eta(1)%init( left=-3.2d0 ,right=3.2d0 ,Nbins=100 )
call h_eta(2)%init( left=-3.2d0 ,right=3.2d0 ,Nbins=100 )

sumW0 = 0
sumW1 = 0
sumW2 = 0

!=======================================================================
do ;call read_event ;if(exitLoop)exit ! ====== DO NOT TOUCH THIS LINE ==
!=======================================================================
!== BEGIN EVENT ========================================================
!=======================================================================

totalWeight = 0

eta(1) = rapidity(pFnst(0:3,1))
eta(2) = rapidity(pFnst(0:3,2))

call sort_big2small( rapOrdrd ,eta(1:2) )

pT(1) = pTrans(pFnst(0:3,rapOrdrd(1)))
if (pT(1).lt.100d0) goto 999

pT(2) = pTrans(pFnst(0:3,rapOrdrd(2)))
flvr(1) = lhefFnst(rapOrdrd(1),procNr)
flvr(2) = lhefFnst(rapOrdrd(2),procNr)

call fragment( flvr(1) ,eta(1),pT(1) ,ww(1),xx(1) )
call fragment( flvr(2) ,eta(2),pT(2) ,ww(2),xx(2) )
pT(1) = xx(1)*pT(1)
pT(2) = xx(2)*pT(2)
totalWeight = eventWeight*ww(1)*ww(2)

999 continue

! histograms must include zero-weight events
call h_pT(1)%collect( pT(1) ,totalWeight )
call h_pT(2)%collect( pT(2) ,totalWeight )
call h_eta(1)%collect( eta(1) ,totalWeight )
call h_eta(2)%collect( eta(2) ,totalWeight )

sumW0 = sumW0+1 ! count number of events
sumW1 = sumW1+totalWeight    ! update sum of weights
sumW2 = sumW2+totalWeight**2 ! update sum of weights**2

!=======================================================================
!== END EVENT ==========================================================
!=======================================================================
enddo ;call close_file ! ===================== DO NOT TOUCH THIS LINE ==
!=======================================================================
!== WRITE HISTOGRAM FILES ==============================================
!=======================================================================

call h_pT(1)%write('pT1.hst')
call h_pT(2)%write('pT2.hst')
call h_eta(1)%write('eta1.hst')
call h_eta(2)%write('eta2.hst')

xsection = sumW1/sumW2
errEst = (sumW2/sumW0-average**2)/(sumW0-1)
write(*,'(a32,e16.8,a4,e16.8)') &
  'Total cross section in nanobarn:',xsection,' +/-',errEst

!=======================================================================
contains ! =================================== DO NOT TOUCH THIS LINE ==
!=======================================================================
!== DEFINE SUBROUTINES =================================================
!=======================================================================

  subroutine fragment( flavor,eta,pT ,ww,xw )
  integer        ,intent( in) :: flavor 
  real(kind(1d0)),intent( in) :: eta,pT
  real(kind(1d0)),intent(out) :: ww,xw
  real(kind(1d0)),parameter :: sqrt2pi=2.5066282746310005024157652848110d0
  real(kind(1d0)),parameter :: conv=5.06952d0
  real(kind(1d0)),parameter :: Lpar=5d0,Kpar=2d0,Qs=0.1d0
  real(kind(1d0)),parameter :: conv3=conv**3
  real(kind(1d0)),parameter :: Lconv2=(Lpar*conv)**2/2
  real(kind(1d0)),parameter :: Lconv3=(Lpar*conv)**3/2
  real(kind(1d0)),parameter :: epsmax=143d0
  real(kind(1d0)),parameter :: a1=2108.05d0,b1=3.66935d0
  real(kind(1d0)),parameter :: a2=486.368d0,b2=1.19377d0
  real(kind(1d0)),parameter :: norm=(a1*b1-a2*b2)*sqrt2pi
  real(kind(1d0)) :: eps,qh,rr,wc,xx,probCont,prob0,jac
  integer :: flvr
!
  eps = epsmax*( a1*exp(-(eta/b1)**2/2) - a2*exp(-(eta/b2)**2/2) )/norm
  qh = 2*Kpar*(eps/conv3)**0.75d0
  rr = qh*Lconv3 ! parameter r - input for the probability
  wc = qh*Lconv2 ! omega_c parameter
  jac = log((pT+Qs)/Qs) ! for the random number mapping
  xx = Qs/wc*( exp(jac*ranfun()) - 1 ) ! energy loss fraction
!
  if (flavor.eq.21) then ! 21=gluon in LHEF encoding
    flvr = 2 ! gluon
  else
    flvr = 1 ! quark
  endif
  call swqmult(flvr,rr,xx,probCont,prob0)
!
  ww = 1
  xw = 1
  if (probCont.ge.0d0) then
    if (ranfun().ge.prob0) then
      xw = (pT-xx*wc)/pT
      ww = (xx+Qs/wc)*jac*max(probCont,0d0)/xw
    endif
  endif
  end subroutine

!=======================================================================
end program ! ================================ DO NOT TOUCH THIS LINE ==
!=======================================================================

!# The 4 collumns in a histogram file represent
!#
!#   left-bin-border  right-bin-border  value  statistical-error-estimate
!#
!# Histogram files are best manipulated with standard unix commands.
!# For example, to put pT distributions without error estimates in one file:
!
!$ paste pT1.hst pT2.hst pT3.hst | awk '{print $1" "$2" "$3" "$7" "$11}' > result.hst
!
!# To get the ratio of pT1 and pT2:
!
!$ paste pT1.hst pT2.hst | awk '{print $1" "$2" "$3" "$7}' \
!$ | awk '{if($4+0!=0) print $1" "$2" "$3/$4; else print $1" "$2" ""0"}' \
!$ > result.hst
!
