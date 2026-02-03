!
! The following obviously-named variables are available:
!
! A size parameter
!   integer,parameter :: Nsize
!
! The number of final states, the number of processes, the number of
! active light flavors:
!   integer :: Nfinst,Nproc,Nflavor
!
! The center-off-mass energy, the positve-rapidity beam energy, and the
! negative-rapidity energy
!   real(kind(1d0)) :: Ecm,EposRap,EnegRap
!
! The process number of the event under consideration
!   integer :: procNr
!
! In the following "A" refers to the positive-rapidity initial state,
! and "B" refers to the negative-rapidity initial state, while "F" refers
! to the final states.
!   integer :: flavorA,colorA,anticA,helicityA
!   integer :: flavorB,colorB,anticB,helicityB
!   integer :: flavorF(Nsize),colorF(Nsize),anticF(Nsize),helicityF(Nsize)
!
! Momentum components are:  0=>E, 1=>px, 2=>py, 3=>pz, 4=>E^2-px^2-py^2-pz^2
! Initial-state momenta have negative energy.
!   real(kind(1d0)) :: pInstA(0:4) ,pInstB(0:4) ,pFinst(0:4,Nsize)
!
! In the following, kT and scale have mass-dimension, not mass-squared
!   real(kind(1d0)) :: pdfB,xB,kTB,scaleB ,pdfA,xA,kTA,scaleA
!
! Furthermore, there is access to
!   real(kind(1d0)) :: eventWeight,matrixElement,partlumi,alphaS,renScale
!
! eventWeight is the only variable that can be changed in this source file.
!
! The following functions are available from the module katie_histogramtools:
!
! real(kind(1d0) function mass(p(0:3))            ! mass = sqrt(abs(square))
! real(kind(1d0) function pTrans(p(0:3))          ! transverse momentum
! real(kind(1d0) function ETrans(p(0:3))          ! transverse mass
! real(kind(1d0) function rapidity(p(0:3))        ! rapidity
! real(kind(1d0) function pseudoRap(p(0:3))       ! pseudo rapidity
! real(kind(1d0) function phi(p(0:3))             ! azimuthal angle
! real(kind(1d0) function theta(p(0:3))           ! polar angle
! real(kind(1d0) function angle(p(0:3),q(0:3))    ! angle between p and q
! real(kind(1d0) function deltaPhi(p(0:3),q(0:3)) ! difference between azim.angles
! real(kind(1d0) function deltaR(p(0:3),q(0:3))   ! delta R
! integer        :: unity(i)   ! the number i 
! character(1)   :: numeral(i) ! the number i as a character
! real(kind(1d0) :: r1PI       ! the number pi
! real(kind(1d0) :: r2PI       ! the number 2*pi
! subroutine sort_big2small(l,a,n) ! sort a(1:n) from big to small or small to big
! subroutine sort_small2big(l,a,n) ! the array l contains the associated permutation
!
! type(breit_type) :: obj ! Lorentz transformation to a Breit frame
!   call obj%init(q(0:3)) ! prepares transformation to Breit frame of q
!   pBreit(0:3) = obj%act(p(0:3)) ! pBreit is the transformed version of p
!
! type(histo_1d_type) :: hst ! 1-dimensional histogram
! type(histo_2d_type) :: hst ! 2-dimensional histogram
!
! Initialize 1-dim histogram with:
!   call hst%init( left=0d0 ,right=2d0 ,Nbins=4 )
! or alternatively:
!   call hst%init( [0.0,0.5 ,0.5,1.0 ,1.0,1.5 ,1.5,2.0] )
! Data collection:
!   call hst%collect( value ,weight )
! Write histogram file:
!   call hst%write('someFileName')
!
! Initialize 2-dim histogram with:
!   call hst%init( xLow=0d0 ,xUpp=2d0 ,NxBins=4 ,yLow=-1d0 ,yUpp=1d0 ,Nybins=2 )
! or alternatively:
!   call hst%init( [0.0,0.5 ,0.5,1.0 ,1.0,1.5 ,1.5,2.0] ,[-1.0,0.0 ,0.0,1.0] )
! 

program create_eventfile !============================ DO NOT TOUCH THIS LINE ==
use katie_eventfile !================================= DO NOT TOUCH THIS LINE ==
use katie_histogramtools !============================ DO NOT TOUCH THIS LINE ==
!== USE OTHER MODULES BELOW THIS LINE ==========================================

implicit none !================= REMOVE THIS LINE AND SUFFER THE CONSEQUENCES ==
!== DECLARE NEW VARIABLES BELOW THIS LINE ======================================

type(histo_1d_type) :: pTa,pTb,pTc,pTd,pTe,pTf
real(kind(1d0)) :: Qsquare,Qvirtual(0:3),pTav
type(breit_type) :: breit


!== CHANGE UNIT NUMBERS IN THE LINE BELOW, IF YOU WISH =========================
call initialize( rawFile_unit=21 ,eventFile_unit=22 ) !=========================
!== INITIALIZE YOUR ROUTINES BELOW THIS LINE ===================================

call pTa%init([5.0,7.0 ,7.0,11.0 ,11.0,18.0 ,18.0,30.0 ,30.0,50.0])
call pTb%init([5.0,7.0 ,7.0,11.0 ,11.0,18.0 ,18.0,30.0 ,30.0,50.0])
call pTc%init([5.0,7.0 ,7.0,11.0 ,11.0,18.0 ,18.0,30.0 ,30.0,50.0])
call pTd%init([5.0,7.0 ,7.0,11.0 ,11.0,18.0 ,18.0,30.0 ,30.0,50.0])
call pTe%init([5.0,7.0 ,7.0,11.0 ,11.0,18.0 ,18.0,30.0 ,30.0,50.0])
call pTf%init([5.0,7.0 ,7.0,11.0 ,11.0,18.0 ,18.0,30.0 ,30.0,50.0])


do ;if (exitLoop>0) exit ;call read_event !=========== DO NOT TOUCH THIS LINE ==
!== COLLECT DATA AND ALTER EVENTWEIGHT BELOW THIS LINE =========================
!== For example: eventWeight = eventWeight/pdfA*my_pdf(flavorA,xA,kTA,scaleA)

Qvirtual(0) = EnegRap-pFinst(0,3)
Qvirtual(1) =        -pFinst(1,3)
Qvirtual(2) =        -pFinst(2,3)
Qvirtual(3) =-EnegRap-pFinst(3,3)
Qsquare = 2*EnegRap*(pFinst(0,3)+pFinst(3,3))
call breit%init(Qvirtual)

pTav = ( pTrans(breit%act(pFinst(0:3,1))) &
        +pTrans(breit%act(pFinst(0:3,2))) )/2

if     ( 150d0.le.Qsquare.and.Qsquare.lt.  200d0) then ;call pTa%collect( pTav ,eventWeight )
elseif ( 200d0.le.Qsquare.and.Qsquare.lt.  270d0) then ;call pTb%collect( pTav ,eventWeight )
elseif ( 270d0.le.Qsquare.and.Qsquare.lt.  400d0) then ;call pTc%collect( pTav ,eventWeight )
elseif ( 400d0.le.Qsquare.and.Qsquare.lt.  700d0) then ;call pTd%collect( pTav ,eventWeight )
elseif ( 700d0.le.Qsquare.and.Qsquare.lt. 5000d0) then ;call pTe%collect( pTav ,eventWeight )
elseif (5000d0.le.Qsquare.and.Qsquare.lt.15000d0) then ;call pTf%collect( pTav ,eventWeight )
endif


call write_event ;enddo !============================= DO NOT TOUCH THIS LINE ==
!== DO STUFF AFTER CREATION OF THE EVENT FILE BELOW THIS LINE ==================

call pTa%write( 'pTa.hst' )
call pTb%write( 'pTb.hst' )
call pTc%write( 'pTc.hst' )
call pTd%write( 'pTd.hst' )
call pTe%write( 'pTe.hst' )
call pTf%write( 'pTf.hst' )


end program !========================================= DO NOT TOUCH THIS LINE ==
