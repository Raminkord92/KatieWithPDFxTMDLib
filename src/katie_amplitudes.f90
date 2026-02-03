!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                              
!! usage of katamp                                                              
!!                                                                              
!! call katamp_add_kinematics( kinID ,Nxtrn ,inStates )                         
!!    output:                                                                   
!!           kinID = integer ID for the kinematics                              
!!    input:                                                                    
!!           Nxtrn = integer number of external particles in the processes      
!!        inStates = integer array of size at least Nxtrn indicating            
!!                   which external particles are initial states,               
!!                     by having a non-zero entry,                              
!!                   and which ones of those are off-shell,                     
!!                     by having the non-zero entry be equal to 2               
!!                                                                              
!!                                                                              
!! call katamp_add_process( kinID ,prcID ,process ,pNonQCD )                    
!!    output:                                                                   
!!           prcID = integer ID for the process                                 
!!    input:                                                                    
!!           kinID = integer ID for the kinematics                              
!!         process = integer array of size at least Nxtrn encoding the process  
!!                   The charges of particles in process should sum up to zero. 
!!         pNonQCD = integer array of size at least 3, containing the powers of 
!!                   non-QCD couplings:                                         
!!                   pNonQCD(3) = power of Higgs-photon couplings               
!!                   pNonQCD(2) = power of Higgs-gluon couplings                
!!                   pNonQCD(1) = power of EW couplings other than              
!!                                the ones mentioned above                      
!!                                                                              
!!                                                                              
!! call katamp_put_momenta( kinID ,momenta ,momSqr )                              
!! call katamp_evaluate( kinID,prcID ,iEval1 ,helicities1 )                     
!! call katamp_evaluate( kinID,prcID ,iEval2 ,helicities2 )                     
!!   ...                                                                        
!! call katamp_evaluate( kinID,prcID ,iEvalN ,helicitiesN )                     
!!                                                                              
!!    Evaluates the color connected amplitudes. For each call to                
!!    katamp_evaluate, amplitudes with the same momenta but different           
!!    helicities are evaluated. The squared amplitude can then be accessed by   
!!    the function katamp_sqr.                                                  
!!    output:                                                                   
!!           iEval = integer id for the evaluated amplitudes                    
!!    input:                                                                    
!!          kinID = integer ID for the kinematics                               
!!          prcID = integer ID for the process                                  
!!        momenta = real array, of shape (0:3,*) with second dimension of size  
!!                  least Nxtrn, containing the values of the momenta of the    
!!                  external  particles in the process. The momenta should sum  
!!                  up to zero.  momenta(0,) is the energy, momenta(1,) is the  
!!                  x-component, etc.                                           
!!         momSqr = real array, of size at least Nxtrn, containing the squared
!!                  momenta (preferably the desired, not the computed, ones.)
!!     helicities = integer array, of size at least Nxtrn, containing the       
!!                  values of the helicities of the external particles in the   
!!                  process. The helicities have one of the values {-1,0,1}, or 
!!                  value WARDTEST which tests the Ward identity.               
!!                                                                              
!! call katamp_sqr( kinID,prcID ,iEval ,rslt )                                  
!!    Returns the squared amplitude summed over color, associated with the      
!!    helicity values of evaluation number iEval.                               
!!                                                                              
!! There are more procedures that may be useful to the user.                    
!! All procedures in this source file starting with "katamp_" are public.       
!!                                                                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module katie_amplitudes
  use avh_prnt
  use avh_ioUnits
  use avh_mathcnst
  use avh_doloops
  use avh_lorentz
  use avh_kskeleton
  use avh_fusiontools
  use avh_buildskel
  use avh_evalskel
  use avh_lagrangian
  use avh_colorflow
  use avh_coloredlgn
  use avh_colormatrix
  use avh_kinematics
  use katie_version
  implicit none

  private
  public :: katamp_add_kinematics ,katamp_add_process ,katamp_prepare_itmd
  public :: katamp_evaluate ,katamp_put_momenta ,katamp_put_internal
  public :: katamp_ampSqr ,katamp_sqr ,katamp_sqr_diag ,katamp_evaluate_test
  public :: katamp_itmd_sqr ,katamp_itmd_sqr_diag 
  public :: katamp_Npartial ,katamp_choose_colcon
  public :: katamp_NoffShell ,katamp_Nxternal ,katamp_offShell ,katamp_offgluon
  public :: katamp_print_itmd ,katamp_set_directions
  public :: vecPerm

  public :: WARDTEST,XDEFINED,POLARANG ! propagated from avh_lorentz

! The following is propagated from avh_lagrangian
! Integer particle identifiers. Anti-particles get a minus-sign.
  public :: eleNeu, eleon,uQuark,dQuark
  public ::  muNeu,  muon,cQuark,sQuark
  public :: tauNeu, tauon,tQuark,bQuark
  public :: Wboson,photon,Zboson,gluon,Higgs
! Parameter and flag setters, eg.:
!   call set_mass_and_width( Zboson ,91.2d0 ,2.5d0 )
!   call set_withWeak( 1 )
  public :: set_mass_and_width ,set_complex_mass_scheme
  public :: set_withQCD ,set_withQED ,set_withWeak
  public :: set_withHiggs ,set_withHG ,set_withHA
  public :: get_anti

  type :: katamp_type
    type(  dsskeleton_type) :: mthr
    type(daughterlist_type) :: dlst
    integer,allocatable :: colorMatrix(:,:),itmdMatrix(:,:)
    integer :: iEval,Npure,Neval=0,Ngluon,Npair,instFlav(2)
    integer,allocatable :: embed(:),pureCC(:)
    complex(kind(1d0)),allocatable :: amp(:,:)
    real(kind(1d0)) :: const,factor
    logical :: offgluon(2)
  end type

  type :: katamp_list_type
    integer :: Nproc=0,pZero=0,Nxtrn=0,Ninst=0,Noff=0
    integer :: inst(2),off(2)
    type(katamp_type),allocatable :: prc(:)
    type(lorentz_type),allocatable :: mom(:)
    real(kind(1d0)) :: factor,energy(2)
  end type

  integer,save :: Nkin=0
  type(katamp_list_type),allocatable,save :: kinema(:)

  real(kind(1d0)),allocatable,save :: pTmp(:,:)
  real(kind(1d0)),allocatable,save :: sTmp(:)
  integer,save :: size_pTmp=0

  real(kind(1d0)),parameter :: defAux(0:3)=&
    [500.000000000000000000000000000000000d0&
    ,376.765346588741920186294009909033775d0&
    ,244.957675573499386700859759002923965d0&
    ,219.188527955943440239756699484160553d0]

  integer,parameter :: vecPerm(3)=[3,1,2]

  logical,save :: initz=.true.

contains


  subroutine katamp_set_directions( pA ,pB )
  real(kind(1d0)),intent(in) :: pA(0:3),pB(0:3)
  if (initz) then
    initz = .false.
    call set_unit('banner',-1)
    call version
    call fill_lagrangian(.true.)
  endif
  call set_real_direction( 1 ,pA ,vecPerm )
  call set_real_direction( 2 ,pB ,vecPerm )
  end subroutine 


  subroutine katamp_add_kinematics( kinID ,Nxtrn ,inStates )
  integer,intent(out) :: kinID
  integer,intent(in) :: Nxtrn
  integer,intent(in) :: inStates(*) ! expect size at least Nxtrn
  type(katamp_list_type),allocatable :: tmpList(:)
  real(kind(1d0)) :: dir(0:3)
  integer :: ii,nn
!
  if (initz) then
    initz = .false.
    call set_unit('banner',-1)
    call version
    call fill_lagrangian(.true.)
    call set_real_direction( 1 ,[-1d0,0d0,0d0,-1d0] ,vecPerm )
    call set_real_direction( 2 ,[-1d0,0d0,0d0, 1d0] ,vecPerm )
  endif
!
  Nkin = Nkin+1
  if (.not.allocated(kinema)) allocate(kinema(1))
  nn = size(kinema,1)
  if (Nkin.gt.nn) then
    allocate(tmpList(1:nn))
    tmpList(1:nn) = kinema(1:nn)
    deallocate(kinema)
    allocate(kinema(1:Nkin+1))
    kinema(1:nn) = tmpList(1:nn)
    deallocate(tmpList)
  endif
!
  kinID = Nkin
!
  associate( kin=>kinema(kinID) )
!
  kin%Nproc = 0
  kin%Nxtrn = Nxtrn
  kin%Ninst = 0
  kin%inst = 0
  do ii=1,Nxtrn
    if (inStates(ii).ne.0) then
      kin%Ninst = kin%Ninst+1
      kin%inst(kin%Ninst) = ii
    endif
  enddo
  if (kin%Ninst.eq.0) then
    if (errru.ge.0) write(errru,*) 'ERROR in katamp_add_kinematics: ' &
      ,'less than 1 initial-state parton.'
    stop
  endif
  kin%Noff = 0
  kin%off = 0
  do ii=1,Nxtrn
    if (inStates(ii).eq.2) then
      kin%Noff = kin%Noff+1
      if (kin%Noff.gt.2) then
        if (errru.ge.0) write(errru,*) 'ERROR in katamp_add_kinematics: ' &
          ,'more than 2 off-shell partons.'
        stop
      endif
      kin%off(kin%Noff) = ii
    endif
  enddo
! tree for calculating all internal momenta
  call update_momTree(kin%Nxtrn-2)
! maximum momentum label
  kin%pZero = base(kin%Nxtrn+1)-1
! allocate momenta
  if (allocated(kin%mom)) deallocate(kin%mom)
  allocate(kin%mom(0:kin%pZero))
!
  end associate
  end subroutine


  subroutine katamp_add_process( kinID ,prcID ,process ,pNonQCD )
  integer,intent(in ) :: kinID
  integer,intent(out) :: prcID
  integer,intent(in) :: process(*) ! expect size at least kinema(kinID)%Nxtrn
  integer,intent(in) :: pNonQCD(*) ! expects size at least 3
  type(katamp_type),allocatable :: tmpList(:)
  integer :: nn
  associate( kin=>kinema(kinID) )
!
  kin%Nproc = kin%Nproc+1
!
  if (.not.allocated(kin%prc)) allocate(kin%prc(2))
  nn = size(kin%prc,1)
  if (kin%Nproc.gt.nn) then
    allocate(tmpList(1:nn))
    tmpList(1:nn) = kin%prc(1:nn)
    deallocate(kin%prc)
    allocate(kin%prc(1:kin%Nproc*2))
    kin%prc(1:nn) = tmpList(1:nn)
    deallocate(tmpList)
  endif
!
  prcID = kin%Nproc
  call init_process( kinID ,prcID ,process(1:kin%Nxtrn) ,pNonQCD )
  if (allocated(kin%prc(prcID)%amp)) deallocate(kin%prc(prcID)%amp)
  allocate(kin%prc(prcID)%amp(1:kin%prc(prcID)%dlst%Ndaughters,1:1))
!
  end associate
  end subroutine


  subroutine init_process( kinID ,prcID ,process ,pNonQCD )
  integer,intent(in) :: kinID,prcID
  integer,intent(in) :: process(1:kinema(kinID)%Nxtrn)
  integer,intent(in) :: pNonQCD(3)
  type(colorContent_type) :: content
  integer :: maxVertMult,ii,Nembed,i1,i2,jj,per(17),rep(17)
  logical :: incomplete,done(17),isPure
  integer,allocatable :: tmpCC(:)
  logical,allocatable :: remove(:)
!
  associate( kin=>kinema(kinID) &
            ,prc=>kinema(kinID)%prc(prcID) ,Nxtrn=>kinema(kinID)%Nxtrn &
            ,Ndght=>kinema(kinID)%prc(prcID)%dlst%Ndaughters )
!
!  if (messu.gt.0) write(messu,*) 'MESSAGE from katamp: ' &
!    ,'preparing 0 ->',trim(prnt_particle(Nxtrn,process))
!
! set the gauge
  call set_gauge( momDef([rZRO,rZRO,rZRO,rZRO],vecPerm) ,rONE )
!
! determine some parameters
  prc%Ngluon = 0
  do ii=1,kin%Nxtrn
    if (process(ii).eq.gluon) prc%Ngluon = prc%Ngluon+1
  enddo
  prc%instFlav(1:2) = process(kin%inst(1:2))
!
! create embedding
  prc%offgluon = .false.
  Nembed = Nxtrn+kin%Noff
  if(allocated(prc%embed))deallocate(prc%embed);allocate(prc%embed(Nembed))
  prc%const = 1
  prc%embed(1:Nxtrn) = process(1:Nxtrn)
  do ii=1,kin%Noff
    i1 = kin%off(ii)
    prc%const = prc%const/gQCD**2
    i2 = prc%embed(i1)
    prc%offgluon(ii) = (i2.eq.gluon)
    select case (i2)
    case (  gluon) ;prc%embed(i1)=-gStar(ii) ;prc%embed(Nxtrn+ii)= gStar(ii)
    case ( uQuark) ;prc%embed(i1)= uStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case (-uQuark) ;prc%embed(i1)=-uStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case ( dQuark) ;prc%embed(i1)= dStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case (-dQuark) ;prc%embed(i1)=-dStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case ( sQuark) ;prc%embed(i1)= sStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case (-sQuark) ;prc%embed(i1)=-sStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case ( cQuark) ;prc%embed(i1)= cStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case (-cQuark) ;prc%embed(i1)=-cStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case ( bQuark) ;prc%embed(i1)= bStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    case (-bQuark) ;prc%embed(i1)=-bStar(ii) ;prc%embed(Nxtrn+ii)= eikPh(ii)
    end select
  enddo
!
  if (messu.gt.0) write(messu,*) 'MESSAGE from katamp: ' &
    ,'preparing 0 ->',trim(prnt_particle(Nembed,prc%embed))
!
! determine color ordering (replaces subroutine fill_content from avh_colorflow)
  call init_colorflow
  if (allocated(content%Ref)) deallocate(content%Ref)
  allocate(content%Ref(2,Nembed))
  content%Nxtrn = Nembed
  content%Ref = 0
  content%Nadjoint = 0
  done = .false.
  i1 = 0
  i2 = 0
! first, off-shell gluons 
  do ii=1,kin%Noff ;if(done(kin%off(ii)))cycle
    if (prc%embed(kin%off(ii)).eq.-gStar(ii)) then
      i1=i1+1 ;content%Ref(1,  Nxtrn+ii )=i1
      i2=i2+1 ;content%Ref(2,kin%off(ii))=i2
      done(kin%off(ii)) = .true.
    endif
  enddo
! then on-shell initial-state gluons 
  do ii=1,kin%Ninst ;if(done(kin%inst(ii)))cycle
    if (prc%embed(kin%inst(ii)).eq.gluon) then
      i1=i1+1 ;content%Ref(1,kin%inst(ii))=i1
      i2=i2+1 ;content%Ref(2,kin%inst(ii))=i2
      content%Nadjoint = content%Nadjoint+1
      done(kin%inst(ii)) = .true.
    endif
  enddo
! then final-state gluons 
  do ii=1,Nxtrn ;if(done(ii))cycle
    if (prc%embed(ii).eq.gluon) then
      i1=i1+1 ;content%Ref(1,ii)=i1
      i2=i2+1 ;content%Ref(2,ii)=i2
      content%Nadjoint = content%Nadjoint+1
      done(ii) = .true.
    endif
  enddo
! then off-shell (anti-)quarks
  do ii=1,kin%Noff ;if(done(kin%off(ii)))cycle
    if     (is_quark(prc%embed(kin%off(ii)))) then
      i1=i1+1 ;content%Ref(1,kin%off(ii))=i1
      done(kin%off(ii)) = .true.
    elseif (is_antiq(prc%embed(kin%off(ii)))) then
      i2=i2+1 ;content%Ref(2,kin%off(ii))=i2
      done(kin%off(ii)) = .true.
    endif
  enddo
! then on-shell initial-state (anti-)quarks
  do ii=1,kin%Ninst ;if(done(kin%inst(ii)))cycle
    if     (is_quark(prc%embed(kin%inst(ii)))) then
      i1=i1+1 ;content%Ref(1,kin%inst(ii))=i1
      done(kin%inst(ii)) = .true.
    elseif (is_antiq(prc%embed(kin%inst(ii)))) then
      i2=i2+1 ;content%Ref(2,kin%inst(ii))=i2
      done(kin%inst(ii)) = .true.
    endif
  enddo
! then final-state (anti-)quarks
  do ii=1,Nxtrn ;if(done(ii))cycle
    if     (is_quark(prc%embed(ii))) then
      i1=i1+1 ;content%Ref(1,ii)=i1
      done(ii) = .true.
    elseif (is_antiq(prc%embed(ii))) then
      i2=i2+1 ;content%Ref(2,ii)=i2
      done(ii) = .true.
    endif
  enddo
  if (i1.ne.i2) then
    if (errru.ge.0) write(errru,*) 'ERROR in katamp: '&
      ,'different number of colors and anti-colors'
    stop
  endif
  content%Npair = i1
  prc%Npair = content%Npair
!
! construct skeletons
  maxVertMult = 4
  if (withHG) maxVertMult = 5
  call dress_skeleton( incomplete ,prc%mthr ,fullTreeLevel(Nembed,maxVertMult) &
                      ,lagrangianRules ,prc%embed )
!  call print_skeleton( prc%mthr ,lagrangianRules ) !DEBUG
  if (incomplete) call not_possible
  call remove_internal_property( incomplete ,prc%mthr ,eikPh )
  if (incomplete) call not_possible
  call update_skeleton( incomplete ,prc%mthr &
                       ,couplingPowerRules ,xCouplingPower(pNonQCD,Nembed) )
  if (incomplete) call not_possible
  call fill_cc_amp( prc%dlst ,prc%mthr ,content )
!  call print_skeleton( prc%dlst ,prc%mthr ,lagrangianRules ) !DEBUG
!
! Remove permutations with more than one cycle for all-gluon amplitudes.
! This didn't happen in  fill_cc_amp  if there are off-shell gluons.
  if (kin%Noff.gt.0.and.prc%Ngluon.eq.kin%Nxtrn) then
    allocate(remove(Ndght)); remove=.false.
    do ii=1,Ndght
      call permutation( per ,prc%Npair ,prc%dlst%daughter(ii)%label )
      i2 = 1
      i1 = 0
      do ;i1=i1+1
        i2 = per(i2)
        if (i2.eq.1) exit
      enddo
      remove(ii) = (i1.lt.prc%Npair)
    enddo
    call prc%dlst%trim(remove)
    deallocate(remove)
  endif
!  call print_skeleton( prc%dlst ,prc%mthr ,lagrangianRules ) !DEBUG
!
! fill color matrix
  allocate(prc%colorMatrix(Ndght,Ndght))
  call colormatrix_square( prc%Npair ,prc%dlst%daughter(:)%Label ,prc%colorMatrix )
!
! Put the amplitude-evaluation counter to zero.
  prc%iEval = 0
!
! Determine list of pure color connections
  allocate(tmpCC(Ndght))
  prc%Npure = 0
  do ii=1,Ndght
    associate( info=>prc%dlst%daughter(ii)%info )
    isPure = .true.
    do jj=1,Nxtrn
      isPure = isPure.and.((info(1,jj).ne.info(2,jj)).or.(info(1,jj).eq.0))
    enddo
    do jj=1,kin%Noff
      if (prc%offgluon(jj)) then
        isPure = isPure.and.(info(2,kin%off(jj)).ne.info(1,Nxtrn+jj))
      endif
    enddo
    if (isPure) then
      prc%Npure = prc%Npure+1
      tmpCC(prc%Npure) = ii
    endif
    end associate
  enddo
  if(allocated(prc%pureCC))deallocate(prc%pureCC);allocate(prc%pureCC(prc%Npure))
  prc%pureCC(1:prc%Npure) = tmpCC(1:prc%Npure)
  deallocate(tmpCC)
!
  end associate
! 
  contains
! 
    subroutine not_possible
    write(*,*) 'ERROR in katamp_add_process: process not possible.'
    stop
    end subroutine
! 
  end subroutine


  subroutine katamp_prepare_itmd( kinID ,prcID ,option )
! Only the following values of  option  are not ignored:
!   option=2 => leading color
  integer,intent(in) :: kinID,prcID ,option
  integer,allocatable :: posLab(:)
  integer :: ii,jj,rUnit
  character(3) :: otherInst,aNpair,aNgluon
  character(254),parameter :: fileName =&
![itmdColorTableFile
!]itmdColorTableFile
!
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
!
  if (kin%Noff.ne.1) call not_possible
  if (.not.prc%offgluon(1)) call not_possible
  if (prc%Npair.gt.6) call not_possible
  if (prc%Npair.eq.6.and.prc%Ngluon.ne.kin%Nxtrn) call not_possible
!
  write(aNpair ,'(i3)') prc%Npair  ;aNpair =adjustl(aNpair ) 
  write(aNgluon,'(i3)') prc%Ngluon ;aNgluon=adjustl(aNgluon) 
  ii=1 ;if(kin%off(ii).eq.kin%inst(ii)) ii=2
  if     (prc%instFlav(ii).eq.gluon)  then ;otherInst='0  '
  elseif (is_quark(prc%instFlav(ii))) then ;otherInst='1  '
  elseif (is_antiq(prc%instFlav(ii))) then ;otherInst='-1 '
  else                                     ;otherInst='2  '
  endif
!
  if (allocated(prc%itmdMatrix))deallocate(prc%itmdMatrix)
  allocate(prc%itmdMatrix(prc%dlst%Ndaughters,prc%dlst%Ndaughters))
!
  allocate(posLab(igamm(prc%Npair+1)));posLab=0
  do ii=1,prc%dlst%Ndaughters
    posLab(prc%dlst%daughter(ii)%label) = ii
  enddo
!
  prc%const = prc%const*(Ncolor(2)-1)
!
  prc%itmdMatrix  = 0
  prc%colorMatrix = 0
!
  open(newunit=rUnit,file=trim(fileName),status='old')
  call print_tablehead
  call read_matrices( &
          trim('BEGIN '//trim(aNpair)//' '//trim(aNgluon)//' '//trim(otherInst)) &
         ,'END' )
  close(rUnit)
!
  do ii=1,prc%dlst%Ndaughters
  do jj=1,prc%dlst%Ndaughters
    if (prc%itmdMatrix(ii,jj).eq.0) then
      prc%colorMatrix(ii,jj) = 0
      prc%itmdMatrix( ii,jj) = 1
    endif
  enddo
  enddo
!
  deallocate(posLab)
  end associate
! 
  contains
!
    subroutine read_matrices( bgnTag ,endTag )
    character(*),intent(in) :: bgnTag,endTag
    character(144) :: line
    integer :: ii,jj,powerN,tmdLabel,pi,pj
    do
      read(rUnit,'(A)',iostat=ii) line
      if (ii.ne.0) then
        write(*,*) 'ERROR in katie_amplitudes: ITMD matrix not available.'
        stop
      endif
      if (line(1:len(bgnTag)).eq.bgnTag) exit
    enddo
    do
      read(rUnit,'(A)') line
      if (line(1:len(endTag)).eq.endTag) exit
      read(line,*) ii,jj,tmdLabel,powerN
      pi = posLab(ii)
      pj = posLab(jj)
      if (pi.ne.0.and.pj.ne.0) then
        kinema(kinID)%prc(prcID)%itmdMatrix( pi,pj) = tmdLabel
        kinema(kinID)%prc(prcID)%colorMatrix(pi,pj) = Ncolor(powerN)
        if (option.eq.2.and.powerN.lt.kinema(kinID)%prc(prcID)%Npair-2) &
          kinema(kinID)%prc(prcID)%colorMatrix(pi,pj) = 0
      endif
    enddo
    end subroutine
! 
    subroutine print_tablehead
    character(144) :: line
    logical,save :: done=.false.
    integer :: ii
    if (done) return ;done = .true.
    do
      read(rUnit,'(A)',iostat=ii) line
      if (ii.ne.0) then
        write(*,*) 'WARNING in katie_amplitudes: ITMD matrix headless.'
        return
      endif
      if (line(1:11).eq.'BEGINheader') exit
    enddo
    do
      read(rUnit,'(A)') line
      if (line(1:9).eq.'ENDheader') exit
      write(*,*) trim(line)
    enddo
    end subroutine
! 
    subroutine not_possible
    write(*,*) 'ERROR in katie_amplitudes: ITMD process not available.'
    stop
    end subroutine
! 
  end subroutine


  subroutine katamp_print_itmd( kinID ,prcID ,wUnit )
  integer,intent(in) :: kinID,prcID,wUnit
  integer :: ii,jj,powerN
  character(256) :: line
  character(3) :: Fn(10)=['qg1','qg2','qg3','gg1','gg2','gg3','gg4','gg5','gg6','gg7']
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
  write(wUnit,'(A)') 'Partial amplitude permutation labels:'
  write(wUnit,'(99i6)') prc%dlst%daughter(1:prc%dlst%Ndaughters)%Label
  write(wUnit,'(A)') 'ITMD color matrix (up to an overall factor) as "PowerOfNc,TMD":'
  do ii=1,prc%dlst%Ndaughters
    line = ''
    do jj=1,prc%dlst%Ndaughters
      if (prc%colorMatrix(ii,jj).eq.0) then
        line = trim(line)//'     0'
      else
        powerN = 0
        do ;if (prc%colorMatrix(ii,jj)/Ncolor(powerN).eq.1) exit
          powerN = powerN+1
        enddo
        line = trim(line)//prnt(powerN,1) &
                         //','//Fn(prc%itmdMatrix(ii,jj))
      endif
    enddo
    write(wUnit,'(A)') trim(line)
  enddo
  end associate
  end subroutine


  function katamp_Npartial( kinID ,prcID ) result(rslt)
  integer,intent(in) :: kinID ,prcID
  integer :: rslt
  rslt = kinema(kinID)%prc(prcID)%dlst%Ndaughters
  end function

  function katamp_Nxternal( kinID ) result(rslt)
  integer,intent(in) :: kinID
  integer :: rslt
  rslt = kinema(kinID)%Nxtrn
  end function

  function katamp_NoffShell( kinID ) result(rslt)
  integer,intent(in) :: kinID
  integer :: rslt
  rslt = kinema(kinID)%Noff
  end function

  function katamp_offShell( kinID ) result(rslt)
  integer,intent(in) :: kinID
  integer :: rslt(2)
  rslt = kinema(kinID)%off
  end function

  function katamp_offgluon( kinID ,prcID ,ii ) result(rslt)
  integer,intent(in) :: kinID,prcID ,ii
  integer :: rslt
  rslt = 0
  if (0.lt.ii.and.ii.le.kinema(kinID)%Noff) then
    if (kinema(kinID)%prc(prcID)%offgluon(ii)) rslt = 1 
  endif
  end function


  subroutine katamp_put_momenta( kinID ,momenta ,momSqr ) !NEW
!OLD  subroutine katamp_put_momenta( kinID ,momenta )
! Put external momenta for amplitude evaluation.
! All internal momneta are computed here.
  integer,intent(in) :: kinID
  real(kind(1d0)),intent(in) :: momenta(0:3,*),momSqr(*) !NEW
!OLD  real(kind(1d0)),intent(in) :: momenta(0:3,*)
  integer :: ii,jj,i1,i2
!
  associate( kin=>kinema(kinID) ,Nxtrn=>kinema(kinID)%Nxtrn )
!
  if (kin%pZero.gt.size_pTmp) then
    size_pTmp = kin%pZero
    if (allocated(pTmp)) deallocate(pTmp)
    if (allocated(sTmp)) deallocate(sTmp)
    allocate(pTmp(0:3,0:size_pTmp))
    allocate(sTmp(    0:size_pTmp))
  endif
!
  pTmp(0:3,0) = 0
  pTmp(0:3,kin%pZero) = 0
  sTmp(0) = 0
  sTmp(kin%pZero) = 0
!
  do ii=1,Nxtrn
    jj = base(ii)
    pTmp(0:3,jj) = momenta(0:3,ii)
    sTmp(    jj) = momSqr(ii) !NEW
!OLD    sTmp(    jj) = ( pTmp(0,jj)+pTmp(3,jj) )*( pTmp(0,jj)-pTmp(3,jj) ) &
!OLD                 - pTmp(1,jj)*pTmp(1,jj) - pTmp(2,jj)*pTmp(2,jj)
    pTmp(0:3,kin%pZero-jj) =-pTmp(0:3,jj)
    sTmp(    kin%pZero-jj) = sTmp(    jj)
  enddo
  do ii=1,momTree(Nxtrn)%N
    i1 = momTree(Nxtrn)%p1(ii)
    i2 = momTree(Nxtrn)%p2(ii)
    jj = i1+i2
    pTmp(0:3,jj) = pTmp(0:3,i1)+pTmp(0:3,i2)
    sTmp(    jj) = ( pTmp(0,jj)+pTmp(3,jj) )*( pTmp(0,jj)-pTmp(3,jj) ) &
                 - pTmp(1,jj)*pTmp(1,jj) - pTmp(2,jj)*pTmp(2,jj)
    pTmp(0:3,kin%pZero-jj) =-pTmp(0:3,jj)
    sTmp(    kin%pZero-jj) = sTmp(    jj)
  enddo
!
  do ii=0,kin%pZero
    kin%mom(ii) = momDef( pTmp(0:3,ii) ,sTmp(ii) ,vecPerm )
  enddo
!
  kin%factor = 1
  do ii=1,kin%Noff
    kin%factor = kin%factor*abs(sTmp(base(kin%off(ii))))
    kin%energy(ii) = abs(pTmp(0,base(kin%off(ii))))
  enddo
!
  end associate
!
  end subroutine


  subroutine katamp_put_internal( kinID ,psPoint )
! Alternative to put_momenta: 
!   put all, external and internal, momenta for amplitude evaluation. 
! Takes psPoint_type as input.
  integer,intent(in) :: kinID
  type(psPoint_type),intent(in) :: psPoint
  integer :: ii
  associate( kin=>kinema(kinID) )
  do ii=0,kin%pZero
    kin%mom(ii) = momDef( psPoint%p(ii)%E,psPoint%p(ii)%V,psPoint%p(ii)%S ,vecPerm )
  enddo
  kin%factor = 1
  do ii=1,kin%Noff
    kin%factor = kin%factor*abs(psPoint%p(base(kin%off(ii)))%S)
    kin%energy(ii) = abs(psPoint%p(base(kin%off(ii)))%E)
  enddo
  end associate
  end subroutine
 
  
  subroutine katamp_evaluate( kinID ,prcID ,iEval ,helicities )
  integer,intent(in) :: kinID,prcID
  integer,intent(out) :: iEval
  integer,intent(in) :: helicities(*)
  type(xternalArgs_type) :: xArg
  type(lorentz_type) :: xCurrent(17)
  integer :: ii,sz1,sz2
  complex(kind(1d0)),allocatable :: tmp(:,:)
!
  associate( kin=>kinema(kinID) &
            ,prc=>kinema(kinID)%prc(prcID) ,Nxtrn=>kinema(kinID)%Nxtrn )
!
  prc%iEval = prc%iEval+1
!
  if (prc%iEval.gt.prc%Neval) then
    prc%Neval = prc%iEval
    sz1 = size(prc%amp,1)
    sz2 = size(prc%amp,2)
    allocate(tmp(1:sz1,1:sz2))
    tmp(1:sz1,1:sz2) = prc%amp(1:sz1,1:sz2)
    deallocate(prc%amp)
    allocate(prc%amp(1:sz1,1:prc%Neval))
    prc%amp(1:sz1,1:sz2) = tmp(1:sz1,1:sz2)
    deallocate(tmp)
  endif
!
  iEval = prc%iEval
!
  do ii=1,Nxtrn
    xArg%identity = prc%embed(ii)
    xArg%helicity = helicities(ii)
!    if (ii.eq.3) xArg%helicity = WARDTEST !DEBUG
    xArg%momentum = kin%mom(base(ii))
    xArg%auxVec = momDef( defAux(0:3) ,vecPerm )
    xCurrent(ii) = lagrangianRules%Particle(xArg%identity)%Xternal(xArg)
  enddo
  do ii=1,kin%Noff
    xArg%identity = prc%embed(Nxtrn+ii)
    xArg%helicity =-helicities(kin%off(ii))
    xArg%momentum = kin%mom(0) !kin%mom(mod(base(Nxtrn+ii),kin%pZero+1))
    xArg%auxVec = momDef( defAux(0:3) ,vecPerm )
    xCurrent(Nxtrn+ii) = lagrangianRules%Particle(xArg%identity)%Xternal(xArg)
  enddo
  call evaluate( prc%amp(:,iEval) ,prc%mthr ,lagrangianRules ,prc%dlst &
                ,kin%mom ,kin%pZero ,xCurrent )
  prc%factor = prc%const
  do ii=1,kin%Noff
    prc%factor = prc%factor*kin%energy(ii)
    if (prc%offgluon(ii)) prc%factor = prc%factor*kin%energy(ii)
  enddo
  end associate
  end subroutine


  subroutine katamp_evaluate_test( kinID ,prcID ,iEval ,helicities ,auxMom ,xTermin )
  integer,intent(in) :: kinID,prcID
  integer,intent(out) :: iEval
  integer,intent(in) :: helicities(*)
  complex(kind(1d0)),intent(in) :: auxMom(0:3,*)
  type(lorentz_type),intent(in) :: xTermin(*)
  type(xternalArgs_type) :: xArg
  type(lorentz_type) :: xCurrent(17)
  integer :: ii,sz1,sz2
  complex(kind(1d0)),allocatable :: tmp(:,:)
!
  associate( kin=>kinema(kinID) &
            ,prc=>kinema(kinID)%prc(prcID) ,Nxtrn=>kinema(kinID)%Nxtrn )
!
  prc%iEval = prc%iEval+1
!
  if (prc%iEval.gt.prc%Neval) then
    prc%Neval = prc%iEval
    sz1 = size(prc%amp,1)
    sz2 = size(prc%amp,2)
    allocate(tmp(1:sz1,1:sz2))
    tmp(1:sz1,1:sz2) = prc%amp(1:sz1,1:sz2)
    deallocate(prc%amp)
    allocate(prc%amp(1:sz1,1:prc%Neval))
    prc%amp(1:sz1,1:sz2) = tmp(1:sz1,1:sz2)
    deallocate(tmp)
  endif
!
  iEval = prc%iEval
!
  do ii=1,Nxtrn
    xArg%identity = prc%embed(ii)
    if     (helicities(ii).eq.-USEAUXMOM) then
      xArg%helicity = 1
    elseif (helicities(ii).eq. USEAUXMOM) then
      xArg%helicity = USEAUXMOM
    else
      xArg%helicity = helicities(ii)
    endif
!    if (ii.eq.1) xArg%helicity = WARDTEST !DEBUG
    xArg%momentum = kin%mom(base(ii))
    xArg%auxVec = momDef(auxMom(0:3,ii) ,vecPerm)
    if (helicities(ii).eq.XDEFINED) then
      xCurrent(ii) = xTermin(ii)
    else
      xCurrent(ii) = lagrangianRules%Particle(xArg%identity)%Xternal(xArg)
    endif
  enddo
  do ii=1,kin%Noff
    xArg%identity = prc%embed(Nxtrn+ii)
    if     (helicities(kin%off(ii)).eq.-USEAUXMOM) then
      xArg%helicity =-USEAUXMOM
    elseif (helicities(kin%off(ii)).eq. USEAUXMOM) then
      xArg%helicity =-1
    else
      xArg%helicity =-helicities(kin%off(ii))
    endif
    xArg%momentum = kin%mom(0) !kin%mom(mod(base(Nxtrn+ii),kin%pZero+1))
    xArg%auxVec = momDef(auxMom(0:3,kin%off(ii)),vecPerm)
    xCurrent(Nxtrn+ii) = lagrangianRules%Particle(xArg%identity)%Xternal(xArg)
  enddo
  call evaluate( prc%amp(:,iEval) ,prc%mthr ,lagrangianRules ,prc%dlst &
                ,kin%mom ,kin%pZero ,xCurrent )
  prc%factor = prc%const
  do ii=1,kin%Noff
    prc%factor = prc%factor*kin%energy(ii)
    if (prc%offgluon(ii)) prc%factor = prc%factor*kin%energy(ii)
  enddo
  end associate
  end subroutine


  subroutine katamp_choose_colcon( kinID ,prcID ,diag ,colcon )
  integer        ,intent(in ) :: kinID ,prcID
  real(kind(1d0)),intent(in ) :: diag(*) ! expects size prc%dlst%Ndaughters
  integer        ,intent(out) :: colcon(2,*) ! expects size kin%Nxtrn
  integer :: ii,imax
  associate( dght=>kinema(kinID)%prc(prcID)%dlst%daughter &
            ,kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
!
  imax = prc%pureCC(maxloc(diag(prc%pureCC),1))
!
  do ii=1,kin%Nxtrn
    colcon(1:2,ii) = prc%dlst%daughter(imax)%info(1:2,ii)
  enddo
  do ii=1,kin%Noff
    if (prc%offgluon(ii)) then
      colcon(1,kin%off(ii)) = prc%dlst%daughter(imax)%info(1,kin%Nxtrn+ii)
    endif
  enddo
  end associate
  end subroutine


  subroutine katamp_ampSqr( kinID ,prcID ,iEval,jEval ,rslt )
  integer,intent(in) :: kinID ,prcID ,iEval,jEval
  complex(kind(1d0)),intent(out) :: rslt
  complex(kind(1d0)) :: zz
  integer :: jj,ii
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
  rslt = 0
  do jj=1,prc%dlst%Ndaughters
    zz = 0
    do ii=1,prc%dlst%Ndaughters
      zz = zz + prc%amp(ii,iEval)*prc%colorMatrix(ii,jj)
    enddo
    rslt = rslt + zz*conjg(prc%amp(jj,jEval))
  enddo
  rslt = rslt*kin%factor*prc%factor
! Put the amplitude-evaluation counter to zero.
  prc%iEval = 0
  end associate
  end subroutine


  subroutine katamp_sqr( kinID ,prcID ,iEval ,rslt )
  integer        ,intent(in ) :: kinID ,prcID ,iEval
  real(kind(1d0)),intent(out) :: rslt
  complex(kind(1d0)) :: zz
  integer :: jj,ii
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
  rslt = 0
  do jj=1,prc%dlst%Ndaughters
    zz = 0
    do ii=1,prc%dlst%Ndaughters
      zz = zz + prc%amp(ii,iEval)*prc%colorMatrix(ii,jj)
    enddo
    rslt = rslt + zz*conjg(prc%amp(jj,iEval))
  enddo
  rslt = rslt*kin%factor*prc%factor
! Put the amplitude-evaluation counter to zero.
  prc%iEval = 0
  end associate
  end subroutine


  subroutine katamp_sqr_diag( kinID ,prcID ,iEval ,rslt ,diag )
  integer         ,intent(in ) :: kinID ,prcID ,iEval
  real(kind(1d0)) ,intent(out) :: rslt ,diag(*) ! expects size prc%dlst%Ndaughters
  complex(kind(1d0)) :: zz
  integer :: jj,ii
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
  rslt = 0
  do jj=1,prc%dlst%Ndaughters
    diag(jj) = prc%amp(jj,iEval)*conjg(prc%amp(jj,iEval)) !*kin%factor*prc%factor
    zz = 0
    do ii=1,prc%dlst%Ndaughters
      if (ii.eq.jj) cycle
      zz = zz + prc%amp(ii,iEval)*prc%colorMatrix(ii,jj)
    enddo
    rslt = rslt + zz*conjg(prc%amp(jj,iEval)) + diag(jj)*prc%colorMatrix(jj,jj)
  enddo
  rslt = rslt*kin%factor*prc%factor
! Put the amplitude-evaluation counter to zero.
  prc%iEval = 0
  end associate
  end subroutine


  subroutine katamp_itmd_sqr( kinID ,prcID ,iEval ,tmds ,rslt )
  integer        ,intent(in ) :: kinID ,prcID ,iEval
  real(kind(1d0)),intent(in ) :: tmds(*) ! expects size maxval(prc%itmdMatrix)
  real(kind(1d0)),intent(out) :: rslt ! expects size prc%dlst%Ndaughters
  complex(kind(1d0)) :: zz
  integer :: jj,ii
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
  rslt = 0
  do jj=1,prc%dlst%Ndaughters
    zz = 0
    do ii=1,prc%dlst%Ndaughters
      zz = zz + prc%amp(ii,iEval)*prc%colorMatrix(ii,jj) &
                                 *tmds(prc%itmdMatrix(ii,jj))
    enddo
    rslt = rslt + zz*conjg(prc%amp(jj,iEval))
  enddo
  rslt = rslt*kin%factor*prc%factor
! Put the amplitude-evaluation counter to zero.
  prc%iEval = 0
!  write(*,'(99e8.1)') abs(prc%amp(:,iEval)) !DEBUG
  end associate
  end subroutine


  subroutine katamp_itmd_sqr_diag( kinID ,prcID ,iEval ,tmds ,rslt ,diag )
  integer        ,intent(in ) :: kinID ,prcID ,iEval
  real(kind(1d0)),intent(in ) :: tmds(*) ! expects size maxval(prc%itmdMatrix)
  real(kind(1d0)),intent(out) :: rslt ,diag(*) ! expects size prc%dlst%Ndaughters
  complex(kind(1d0)) :: zz
  integer :: jj,ii
  associate( kin=>kinema(kinID) ,prc=>kinema(kinID)%prc(prcID) )
  rslt = 0
  do jj=1,prc%dlst%Ndaughters
    diag(jj) = prc%amp(jj,iEval)*conjg(prc%amp(jj,iEval)) !*kin%factor*prc%factor
    zz = 0
    do ii=1,prc%dlst%Ndaughters
      if (ii.eq.jj) cycle
      zz = zz + prc%amp(ii,iEval)*prc%colorMatrix(ii,jj) &
                                 *tmds(prc%itmdMatrix(ii,jj))
    enddo
    rslt = rslt + zz*conjg(prc%amp(jj,iEval)) &
                + diag(jj)*prc%colorMatrix(jj,jj)*tmds(prc%itmdMatrix(jj,jj))
  enddo
  rslt = rslt*kin%factor*prc%factor
! Put the amplitude-evaluation counter to zero.
  prc%iEval = 0
  end associate
  end subroutine


end module


