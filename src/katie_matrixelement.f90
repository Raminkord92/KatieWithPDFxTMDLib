module katie_matrixelement
  use avh_prnt
  use avh_ioUnits
  use avh_kskeleton, only: size_xternal
  use avh_ranXconf
  use avh_lagrangian
  use avh_kinematics
  use katie_amplitudes
  implicit none

  private
  public :: add_kinematics ,matrixelement_type

  type :: matrixelement_type
    integer :: kinID,prcID,Nxtrn,Nprtl,Noff,off(2)
    type(ranXconf_type) :: xHelConf,osHelConf
    real(kind(1d0)) :: const
  contains
    procedure :: init
    procedure :: hSum_1 
    procedure :: hSum_2
    procedure :: hSum_3
    procedure :: hSum_4
    procedure :: hSum_5
    procedure :: hSum_6
    generic :: helicitySummed=>hSum_1,hSum_2,hSum_3,hSum_4,hSum_5,hSum_6
    procedure :: hSampl_1
    procedure :: hSampl_2
    procedure :: hSampl_3
    procedure :: hSampl_4
    generic :: helicitySampled=>hSampl_1,hSampl_2,hSampl_3,hSampl_4
    procedure :: print_itmd
  end type

contains


  subroutine add_kinematics( kinID ,Nxtrn ,inStates )
  integer,intent(out) :: kinID
  integer,intent(in ) :: Nxtrn
  integer,intent(in ) :: inStates(Nxtrn)
  call katamp_add_kinematics( kinID ,Nxtrn ,inStates )
  end subroutine


  subroutine init( obj ,kinID ,process ,pNonQCD ,heltype &
                  ,itmdf ,leadingColor )
!***********************************************************************
!***********************************************************************
  class(matrixelement_type) :: obj
  integer,intent(in) :: kinID
  integer,intent(in) :: process(*) ! expects size Nxtrn
  integer,intent(in) :: pNonQCD(3)
  character(*),intent(in) :: helType
  logical,intent(in),optional :: itmdf,leadingColor
  integer :: ii,xNspinDof(size_xternal),option
!
  obj%kinID = kinID
  call katamp_add_process( obj%kinID ,obj%prcID ,process ,pNonQCD )
!
  obj%Nxtrn = katamp_Nxternal( obj%kinID )
  obj%Nprtl = katamp_Npartial( obj%kinID ,obj%prcID )
  obj%Noff  = katamp_NoffShell( obj%kinID )
  obj%off   = katamp_offShell( obj%kinID )
!
! Correction factor for summing over "helicities"(=spinor-choices) of off-shell gluons.
! The ranXconf_type corrects for unnecessary duplicate calculations regarding
! helicity configurations in general.
  obj%const = 1
  if (katamp_offgluon(obj%kinID,obj%prcID,1).ne.0) obj%const = obj%const/2
  if (katamp_offgluon(obj%kinID,obj%prcID,2).ne.0) obj%const = obj%const/2
!
  xNspinDof(1:obj%Nxtrn) = get_NspinDof(obj%Nxtrn,process)
  select case (helType)
    case ('sum')
      call obj%xHelConf%fill( xNspinDof(1:obj%Nxtrn) ,helVal ,option='sum' )
    case default
! Always explicitly sum over spinor-choices for off-shell partons, also in the
! random-helicity case.
      do ii=1,obj%Noff
        xNspinDof(obj%off(ii)) = 1
      enddo
      call obj%xHelConf%fill( xNspinDof(1:obj%Nxtrn) ,helVal )
      xNspinDof(1:obj%Nxtrn) = 1
      do ii=1,obj%Noff
        xNspinDof(obj%off(ii)) = 2
      enddo
      call obj%osHelConf%fill( xNspinDof(1:obj%Nxtrn) ,helVal )
  end select
!
  if(present(itmdf))then;if(itmdf)then
    option = 0
    if(present(leadingColor))then;if(leadingColor)then
      option = 2
    endif;endif
    call katamp_prepare_itmd( obj%kinID ,obj%prcID ,option )
  endif;endif
!
  end subroutine


  subroutine print_itmd( obj ,wUnit )
  class(matrixelement_type) :: obj
  integer,intent(in) :: wUnit
  call katamp_print_itmd( obj%kinID ,obj%prcID ,wUnit )
  end subroutine


  include 'katie_matrixelement.h90'
  !execute case: itmdf | no
  !execute case: task | optimize
  !execute case: momenta | psPoint

  include 'katie_matrixelement.h90'
  !execute case: itmdf | no
  !execute case: task | mainMC
  !execute case: momenta | psPoint

  include 'katie_matrixelement.h90'
  !execute case: itmdf | yes
  !execute case: task | optimize
  !execute case: momenta | psPoint

  include 'katie_matrixelement.h90'
  !execute case: itmdf | yes
  !execute case: task | mainMC
  !execute case: momenta | psPoint

  include 'katie_matrixelement.h90'
  !execute case: itmdf | no
  !execute case: task | optimize
  !execute case: momenta | array

  include 'katie_matrixelement.h90'
  !execute case: itmdf | yes
  !execute case: task | optimize
  !execute case: momenta | array


!  subroutine helicitySummed( obj ,rslt ,psPoint ,colcon )
!  class(matrixelement_type) :: obj
!  real(kind(1d0))   ,intent(out) :: rslt
!  type(psPoint_type),intent(in ) :: psPoint
!  integer           ,intent(out),optional :: colcon(2,obj%Nxtrn)
!  integer :: helicity(size_xternal),iEval
!  real(kind(1d0)) :: wght,term,diag(720),sumDiag(720)
!  logical :: rtrnColcon
!!
!  rtrnColcon = present(colcon)
!  if (rtrnColcon) sumDiag(1:obj%Nprtl) = 0
!!
!  call katamp_put_internal( obj%kinID ,psPoint )
!  rslt = 0
!  do ;if (obj%xHelConf%exit(helicity(1:obj%Nxtrn),wght)) exit
!    call katamp_evaluate( obj%kinID,obj%prcID ,iEval ,helicity )
!    call katamp_sqr_and_diag( obj%kinID,obj%prcID ,iEval ,term ,diag )
!    call obj%xHelConf%sum( term )
!    rslt = rslt + term*wght
!    if (rtrnColcon) then
!      sumDiag(1:obj%Nprtl) = sumDiag(1:obj%Nprtl) + diag(1:obj%Nprtl)
!    endif
!  enddo
!  rslt = rslt*obj%const
!  if (rtrnColcon) call katamp_choose_colcon( obj%kinID,obj%prcID ,sumDiag ,colcon )
!  end subroutine
!
!
!  subroutine helicitySampled( obj ,rslt ,psPoint ,rho ,colcon ,helicity )
!  class(matrixelement_type) :: obj
!  real(kind(1d0))   ,intent(out) :: rslt
!  type(psPoint_type),intent(in ) :: psPoint
!  real(kind(1d0))   ,intent(in ) :: rho
!  integer  ,intent(out),optional :: colcon(2,obj%Nxtrn),helicity(obj%Nxtrn)
!  integer :: osHelicity(size_xternal),xHelicity(size_xternal),iEval,ii
!  real(kind(1d0)) :: wght,term,diag(720),sumDiag(720)
!  logical :: rtrnColcon
!!
!  call obj%xHelConf%gnrt( xHelicity ,rho )
!!
!  rtrnColcon = present(colcon)
!  if (rtrnColcon) sumDiag(1:obj%Nprtl) = 0
!!
!  call katamp_put_internal( obj%kinID ,psPoint )
!  rslt = 0
!  do ;if (obj%osHelConf%exit(osHelicity(1:obj%Nxtrn),wght)) exit
!    do ii=1,obj%Noff
!      xHelicity(obj%off(ii)) = osHelicity(obj%off(ii))
!    enddo
!    call katamp_evaluate( obj%kinID,obj%prcID ,iEval ,xHelicity )
!    call katamp_sqr_and_diag( obj%kinID,obj%prcID ,iEval ,term ,diag )
!    call obj%osHelConf%sum( term )
!    rslt = rslt + term*wght
!    if (rtrnColcon) sumDiag(1:obj%Nprtl) = sumDiag(1:obj%Nprtl) + diag(1:obj%Nprtl)
!  enddo
!!
!  wght = obj%xHelConf%wght()  
!  rslt = rslt*obj%const
!  call obj%xHelConf%digital( rslt )
!  rslt = rslt*wght
!!
!  if (rtrnColcon) call katamp_choose_colcon( obj%kinID,obj%prcID ,sumDiag ,colcon )
!  if (present(helicity)) helicity(1:obj%Nxtrn) = xHelicity(1:obj%Nxtrn)
!  end subroutine


end module


