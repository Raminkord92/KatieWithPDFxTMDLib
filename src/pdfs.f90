module katie_pdfs
  use avh_iounits
  use avh_lagrangian
  use avh_trivinp
  use avh_mathcnst
  use katie_partlumi ,only: set_Nflavors
  use iso_c_binding ,only: c_ptr,c_null_ptr,c_char,c_null_char,c_int,c_double,c_associated
  implicit none
  private
  public :: pdfs_init,pdfs_init_base,pdfs_init_lhapdf,set_tmdpath,add_tmdpdf,pdfMesh
  public :: pdfvec_lhapdf,pdfvec_grid
  public :: pdfunc_lhapdf,pdfunc_grid
  public :: pdfs_init_tmdlib,pdfvec_tmdlib,pdfunc_tmdlib !{TMDlibVersion!=000!}
  public :: pdfvec_tmdlibset,pdfunc_tmdlibset !{TMDlibVersion=2xx!}
  public :: pdfs_init_pdfxtmd,pdfvec_pdfxtmd,pdfunc_pdfxtmd
  public :: pdfvec_pdfxtmdset,pdfunc_pdfxtmdset
  public :: pdfs_init_pdfxtmd_cpdf,pdfvec_pdfxtmd_cpdf,pdfunc_pdfxtmd_cpdf
  public :: pdfs_init_pdfxtmd_coupling
  public :: alphasFunc

  type(trivinp_type),save,protected :: pdfMesh( -6 -2*13 : 6 -1*13 )
  integer,save :: lgn2pdf(nulParticle:maxParticle)

  integer,parameter :: gl=0,uQ=1,dQ=2,sQ=3,cQ=4,bQ=5,tQ=6
  integer,parameter :: pdf2lha(-6:6)=[-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6]

  integer,save :: tmdFlavor(13),tmdNflavors=0

  type pdfx_state_type
    character(256) :: setLabel = ''
    integer :: member = 0
    type(c_ptr) :: tmdHandle = c_null_ptr
    logical :: isReady = .false.
  end type
  type(pdfx_state_type),save :: pdfx_state(2)
  type pdfx_cpdf_state_type
    character(256) :: setLabel = ''
    integer :: member = 0
    type(c_ptr) :: cpdfHandle = c_null_ptr
    logical :: isReady = .false.
  end type
  type(pdfx_cpdf_state_type),save :: pdfx_cpdf_state(2)
  logical,save :: pdfx_warned = .false.
  type(c_ptr),save :: pdfx_tmd_factory = c_null_ptr
  type(c_ptr),save :: pdfx_cpdf_factory = c_null_ptr
  type(c_ptr),save :: pdfx_coupling_factory = c_null_ptr
  type(c_ptr),save :: pdfx_coupling_handle = c_null_ptr
  character(256),save :: pdfx_coupling_set = ''
  logical,save :: pdfx_use_coupling = .false.

  character(256),save :: path

  interface
    function pdfx_create_coupling_factory() bind(C,name="create_coupling_factory") result(ptr)
      import :: c_ptr
      type(c_ptr) :: ptr
    end function
    subroutine pdfx_delete_coupling_factory(factory) bind(C,name="delete_coupling_factory")
      import :: c_ptr
      type(c_ptr),value :: factory
    end subroutine
    function pdfx_mk_coupling(factory,setName) bind(C,name="mk_coupling") result(ptr)
      import :: c_ptr,c_char
      type(c_ptr),value :: factory
      character(kind=c_char),dimension(*) :: setName
      type(c_ptr) :: ptr
    end function
    function pdfx_alphaQCDMu2_wrapper(cpl,mu2) bind(C,name="alphaQCDMu2_wrapper") result(val)
      import :: c_ptr,c_double
      type(c_ptr),value :: cpl
      real(c_double),value :: mu2
      real(c_double) :: val
    end function
    function pdfx_create_tmd_factory() bind(C,name="create_tmd_factory") result(ptr)
      import :: c_ptr
      type(c_ptr) :: ptr
    end function
    subroutine pdfx_delete_tmd_factory(factory) bind(C,name="delete_tmd_factory")
      import :: c_ptr
      type(c_ptr),value :: factory
    end subroutine
    function pdfx_mk_tmd(factory,setName,setMember) bind(C,name="mk_tmd") result(ptr)
      import :: c_ptr,c_int,c_char
      type(c_ptr),value :: factory
      character(kind=c_char),dimension(*) :: setName
      integer(c_int),value :: setMember
      type(c_ptr) :: ptr
    end function
    function pdfx_tmd_wrapper(tmd,flavor,x,kt2,mu2) bind(C,name="tmd_wrapper") result(val)
      import :: c_ptr,c_int,c_double
      type(c_ptr),value :: tmd
      integer(c_int),value :: flavor
      real(c_double),value :: x,kt2,mu2
      real(c_double) :: val
    end function
    function pdfx_create_cpdf_factory() bind(C,name="create_cpdf_factory") result(ptr)
      import :: c_ptr
      type(c_ptr) :: ptr
    end function
    subroutine pdfx_delete_cpdf_factory(factory) bind(C,name="delete_cpdf_factory")
      import :: c_ptr
      type(c_ptr),value :: factory
    end subroutine
    function pdfx_mk_cpdf(factory,setName,setMember) bind(C,name="mk_cpdf") result(ptr)
      import :: c_ptr,c_int,c_char
      type(c_ptr),value :: factory
      character(kind=c_char),dimension(*) :: setName
      integer(c_int),value :: setMember
      type(c_ptr) :: ptr
    end function
    function pdfx_cpdf_wrapper(cpdf,flavor,x,mu2) bind(C,name="cpdf_wrapper") result(val)
      import :: c_ptr,c_int,c_double
      type(c_ptr),value :: cpdf
      integer(c_int),value :: flavor
      real(c_double),value :: x,mu2
      real(c_double) :: val
    end function
  end interface

!{TMDlibVersion=1,2
  interface pdfs_init_tmdlib
    module procedure init_tmdlib,init_tmdlib_nr
  end interface
!}TMDlibVersion

contains


  subroutine pdfx_ensure_factory()
  if (.not. c_associated(pdfx_tmd_factory)) then
    pdfx_tmd_factory = pdfx_create_tmd_factory()
    if (.not. c_associated(pdfx_tmd_factory)) then
      call pdfx_placeholder_notice('create_tmd_factory failed',1)
    end if
  end if
  end subroutine

  subroutine pdfx_ensure_cpdf_factory()
  if (.not. c_associated(pdfx_cpdf_factory)) then
    pdfx_cpdf_factory = pdfx_create_cpdf_factory()
    if (.not. c_associated(pdfx_cpdf_factory)) then
      call pdfx_placeholder_notice('create_cpdf_factory failed',1)
    end if
  end if
  end subroutine

  subroutine pdfx_ensure_coupling_factory()
  if (.not. c_associated(pdfx_coupling_factory)) then
    pdfx_coupling_factory = pdfx_create_coupling_factory()
    if (.not. c_associated(pdfx_coupling_factory)) then
      call pdfx_placeholder_notice('create_coupling_factory failed',1)
    end if
  end if
  end subroutine

  subroutine pdfx_make_cstring(src,cbuf)
  character(*),intent(in) :: src
  character(kind=c_char),allocatable,intent(out) :: cbuf(:)
  integer :: n,ii
  n = len_trim(src)
  allocate(cbuf(n+1))
  do ii=1,n
    cbuf(ii) = transfer(src(ii:ii),cbuf(ii))
  end do
  cbuf(n+1) = c_null_char
  end subroutine

  integer function pdfx_parton_code(iParton) result(pid)
  integer,intent(in) :: iParton
  select case (iParton)
  case(gl)
    pid = 21
  case(uQ)
    pid = 2
  case(-uQ)
    pid = -2
  case(dQ)
    pid = 1
  case(-dQ)
    pid = -1
  case(sQ)
    pid = 3
  case(-sQ)
    pid = -3
  case(cQ)
    pid = 4
  case(-cQ)
    pid = -4
  case(bQ)
    pid = 5
  case(-bQ)
    pid = -5
  case(tQ)
    pid = 6
  case(-tQ)
    pid = -6
  case default
    pid = 0
  end select
  end function

  integer function pdfx_match_state(label,member,defaultIdx) result(idx)
  character(*),intent(in) :: label
  integer,intent(in) :: member,defaultIdx
  integer :: ii
  idx = 0
  do ii=1,size(pdfx_state)
    if (pdfx_state(ii)%isReady) then
      if (trim(pdfx_state(ii)%setLabel)==trim(label) .and. pdfx_state(ii)%member==member) then
        idx = ii
        return
      end if
    end if
  end do
  if (defaultIdx>0) idx = defaultIdx
  end function

  function pdfx_eval_tmd(idx,xx,kTsq,qq) result(list)
  integer,intent(in) :: idx
  !(realknd2!),intent(in) :: xx,kTsq,qq
  !(realknd2!) :: list(-6:6)
  real(kind(1d0)) :: locX,locKt,locQ,val
  integer :: lf,jj,useIdx
  integer :: pid
  list = 0
  useIdx = idx
  if (.not. pdfx_state(useIdx)%isReady) then
    do jj=1,size(pdfx_state)
      if (pdfx_state(jj)%isReady) then
        useIdx = jj
        exit
      end if
    end do
  end if
  if (.not. pdfx_state(useIdx)%isReady) then
    call pdfx_placeholder_notice('pdfx_eval_tmd init missing',idx)
    return
  end if
  if (.not. c_associated(pdfx_state(useIdx)%tmdHandle)) then
    call pdfx_placeholder_notice('pdfx_eval_tmd handle missing',idx)
    return
  end if
  locX = xx
  if (locX.le.0d0) return
  locKt = max(kTsq,1d-8)
  locQ = max(qq,1d-8)
  do lf=-6,6
    pid = pdfx_parton_code(lf)
    if (pid.ne.0) then
      val = pdfx_tmd_wrapper(pdfx_state(useIdx)%tmdHandle,pid,locX,locKt,locQ)
      ! DEBUG: remove this logging once PDFxTMD validation is done
     ! if (errru.ge.0) write(errru,'(A,1X,I3,1X,3(1X,ES14.6),1X,ES14.6)') &
       ! 'PDFxTMD dbg lf x kt2 mu2 val:',lf,locX,locKt,locQ,val
      list(lf) = val/locX
    end if
  end do
  end function

  function pdfx_eval_cpdf(idx,xx,qq) result(list)
  integer,intent(in) :: idx
  !(realknd2!),intent(in) :: xx,qq
  !(realknd2!) :: list(-6:6)
  real(kind(1d0)) :: locX,locQ,val
  integer :: lf,jj,useIdx
  integer :: pid
  list = 0
  useIdx = idx
  if (.not. pdfx_cpdf_state(useIdx)%isReady) then
    do jj=1,size(pdfx_cpdf_state)
      if (pdfx_cpdf_state(jj)%isReady) then
        useIdx = jj
        exit
      end if
    end do
  end if
  if (.not. pdfx_cpdf_state(useIdx)%isReady) then
    call pdfx_placeholder_notice('pdfx_eval_cpdf init missing',idx)
    return
  end if
  if (.not. c_associated(pdfx_cpdf_state(useIdx)%cpdfHandle)) then
    call pdfx_placeholder_notice('pdfx_eval_cpdf handle missing',idx)
    return
  end if
  locX = xx
  if (locX.le.0d0) return
  locQ = max(qq,1d-8)
  do lf=-6,6
    pid = pdfx_parton_code(lf)
    if (pid.ne.0) then
      val = pdfx_cpdf_wrapper(pdfx_cpdf_state(useIdx)%cpdfHandle,pid,locX,locQ*locQ)
      list(lf) = val/locX
    end if
  end do
  end function


  subroutine pdfs_init( Nflavors ,lhaSet )
  integer,intent(in) :: Nflavors
  character(*),intent(in) :: lhaSet
  call pdfs_init_base( Nflavors )
  call pdfs_init_lhapdf( lhaSet )
  end subroutine

  subroutine pdfs_init_base( Nflavors )
  integer,intent(in) :: Nflavors
  call set_Nflavors( Nflavors )
  lgn2pdf = gl
  lgn2pdf(uQuark)=uQ; lgn2pdf(-uQuark)=-lgn2pdf(uQuark);
  lgn2pdf(dQuark)=dQ; lgn2pdf(-dQuark)=-lgn2pdf(dQuark);
  lgn2pdf(cQuark)=cQ; lgn2pdf(-cQuark)=-lgn2pdf(cQuark);
  lgn2pdf(sQuark)=sQ; lgn2pdf(-sQuark)=-lgn2pdf(sQuark);
  lgn2pdf(tQuark)=tQ; lgn2pdf(-tQuark)=-lgn2pdf(tQuark);
  lgn2pdf(bQuark)=bQ; lgn2pdf(-bQuark)=-lgn2pdf(bQuark);
  tmdNflavors = 0
  end subroutine

  subroutine pdfs_init_lhapdf( lhaSet )
  character(*),intent(in) :: lhaSet
  call InitPDFsetByName(trim(lhaSet)//'.LHgrid')
  end subroutine

  subroutine set_tmdpath( pathVal )
  character(*),intent(in) :: pathVal
  path = pathVal
  end subroutine

  subroutine add_tmdpdf( iBeam ,iParton ,fileName )
  integer,intent(in) :: iBeam ,iParton
  character(*),intent(in) :: fileName
  if (all(tmdFlavor(1:tmdNflavors).ne.lgn2pdf(iParton))) then
    tmdNflavors = tmdNflavors+1
    tmdFlavor(tmdNflavors) = lgn2pdf(iParton)
  endif
  call pdfMesh(lgn2pdf(iParton)+13*iBeam)%read_file(trim(path)//trim(fileName))
  end subroutine

  function alphasFunc(xx) result(rslt)
  !(realknd2!),intent(in) :: xx
  !(realknd2!) :: rslt,alphasPDF
  if (pdfx_use_coupling .and. c_associated(pdfx_coupling_handle)) then
    rslt = pdfx_alphaQCDMu2_wrapper(pdfx_coupling_handle,xx*xx)
  else
    rslt = alphasPDF(xx)
  endif
  end function
  
  function pdfvec_lhapdf( xx ,qq ) result(list)
  intent(in) :: xx,qq
  !(realknd2!) :: xx,qq,list(-6:6)
  list = 0
  call evolvePDF(xx,qq,list)
  list = list(pdf2lha)/xx
  end function

  function pdfunc_lhapdf( iParton ,xx ,qq ) result(rslt)
  intent(in) :: iParton,xx,qq
  integer :: iParton
  !(realknd2!) :: xx,qq,rslt,list(-6:6)
  rslt = 0
  call evolvePDF(xx,qq,list)
  rslt = list(pdf2lha(lgn2pdf(iParton)))/xx
  end function


  function pdfvec_grid( iBeam ,xx ,qq ,kTsq ) result(list)
  intent(in) :: iBeam ,xx,kTsq,qq
  !(realknd2!) :: xx,kTsq,qq,list(-6:6),logx,logkTsq,logQsq,Qsq
  integer :: iBeam,ii
  list = 0
  select case (pdfMesh(tmdFlavor(1)+13*iBeam)%Ndim)
  case default
    logx = log(xx)
    logkTsq = log(kTsq)
    logQsq = log(qq*qq)
    do ii=1,tmdNflavors
      list(tmdFlavor(ii)) = pdfMesh(tmdFlavor(ii)+13*iBeam)%evaluate( logx ,logkTsq ,logQsq )
      list(tmdFlavor(ii)) = list(tmdFlavor(ii))/xx
    enddo
  case (2)
    logx = log(xx)
    logkTsq = log(kTsq)
    do ii=1,tmdNflavors
      list(tmdFlavor(ii)) = pdfMesh(tmdFlavor(ii)+13*iBeam)%evaluate( logx ,logkTsq )
      list(tmdFlavor(ii)) = list(tmdFlavor(ii))/xx
    enddo
  end select
  end function

  function pdfunc_grid( iBeam ,iParton ,xx ,qq ,kTsq ) result(rslt)
  intent(in) :: iBeam,iParton,xx,qq,kTsq
  integer :: iBeam,iParton
  !(realknd2!) :: xx,kTsq,qq,rslt,logx,logkTsq,logQsq,Qsq
  rslt = 0
  select case (pdfMesh(lgn2pdf(iParton)+13*iBeam)%Ndim)
  case default
    logx = log(xx)
    logkTsq = log(kTsq)
    logQsq = log(qq*qq)
    rslt = pdfMesh(lgn2pdf(iParton)+13*iBeam)%evaluate( logx ,logkTsq ,logQsq )
    rslt = rslt/xx
  case (2)
    logx = log(xx)
    logkTsq = log(kTsq)
    rslt = pdfMesh(lgn2pdf(iParton)+13*iBeam)%evaluate( logx ,logkTsq )
    rslt = rslt/xx
  end select
  end function


!{TMDlibVersion=1,2
  subroutine init_tmdlib( tmdSet )
  character(*),intent(in) :: tmdSet
  integer :: TMDnumberPDF
  call TMDinit(TMDnumberPDF(trim(adjustl(tmdSet))//char(0)))
  end subroutine

  subroutine init_tmdlib_nr( tmdSet )
  integer,intent(in) :: tmdSet
  call TMDinit(tmdSet)
  end subroutine

  function pdfvec_tmdlib( xx ,qq ,kTsq ,kf ) result(list)
  intent(in) :: kf,xx,kTsq,qq
  integer :: kf
  !(realknd2!) :: xx,xbar,kTsq,qq,list(-6:6),hh(2)
  list = 0
  xbar = 0
  call TMDpdf( kf ,xx ,xbar ,sqrt(kTsq) ,qq &
      ,list(uQ),list(-uQ) ,list(dQ),list(-dQ) ,list(sQ),list(-sQ) &
      ,list(cQ),list(-cQ) ,list(bQ),list(-bQ) ,list(gl) )
!  if (kf.gt.0) then
!    list(1) = list(1) + list(-1)
!    list(2) = list(2) + list(-2)
!  else
!    hh(1) = list(-1)
!    hh(2) = list(-2)
!    list(-1) = list(1) + list(-1)
!    list(-2) = list(2) + list(-2)
!    list(1:2) = hh(1:2)
!  endif
!  list(-3) = list(3)
!  list(-4) = list(4)
!  list(-5) = list(5)
  list = list/xx
  end function

  function pdfunc_tmdlib( iParton ,xx ,qq ,kTsq ,kf ) result(rslt)
  intent(in) :: iParton,kf,xx,kTsq,qq
  integer :: iParton,kf
  !(realknd2!) :: xx,kTsq,qq,list(-6:6),hh,rslt
  list = pdfvec_tmdlib( xx ,qq ,kTsq ,kf )
  rslt = list(lgn2pdf(iParton))
  end function
!}TMDlibVersion

!{TMDlibVersion=2
  function pdfvec_tmdlibset( xx ,qq ,kTsq ,kf ,iSet ) result(list)
  intent(in) :: iSet,kf,xx,kTsq,qq
  integer :: kf,iSet
  !(realknd2!) :: xx,xbar,kTsq,qq,list(-6:6),hh(2)
  list = 0
  xbar = 0
  call TMDpdfset( iSet ,kf ,xx ,xbar ,sqrt(kTsq) ,qq &
      ,list(uQ),list(-uQ) ,list(dQ),list(-dQ) ,list(sQ),list(-sQ) &
      ,list(cQ),list(-cQ) ,list(bQ),list(-bQ) ,list(gl) )
  list = list/xx
  end function

  function pdfunc_tmdlibset( iParton ,xx ,qq ,kTsq ,kf ,iSet ) result(rslt)
  intent(in) :: iSet,iParton,kf,xx,kTsq,qq
  integer :: iParton,kf,iSet
  !(realknd2!) :: xx,kTsq,qq,list(-6:6),hh,rslt
  list = pdfvec_tmdlibset( xx ,qq ,kTsq ,kf ,iSet )
  rslt = list(lgn2pdf(iParton))
  end function
!}TMDlibVersion

  subroutine pdfs_init_pdfxtmd( setLabel ,setMember ,beamID )
  character(*),intent(in) :: setLabel
  integer,intent(in),optional :: setMember,beamID
  integer :: idx,memberVal
  character(kind=c_char),allocatable :: cset(:)
  idx = 1
  if (present(beamID)) idx = max(1,min(size(pdfx_state),beamID))
  memberVal = 0
  if (present(setMember)) memberVal = setMember
  call pdfx_ensure_factory()
  if (.not. c_associated(pdfx_tmd_factory)) then
    call pdfx_placeholder_notice('PDFxTMD factory missing',idx)
    return
  endif
  call pdfx_make_cstring(setLabel,cset)
  pdfx_state(idx)%tmdHandle = pdfx_mk_tmd(pdfx_tmd_factory,cset,memberVal)
  deallocate(cset)
  if (.not. c_associated(pdfx_state(idx)%tmdHandle)) then
    pdfx_state(idx)%isReady = .false.
    call pdfx_placeholder_notice('mk_tmd returned null',idx)
    return
  endif
  pdfx_state(idx)%setLabel = trim(setLabel)
  pdfx_state(idx)%member = memberVal
  pdfx_state(idx)%isReady = .true.
  end subroutine

  subroutine pdfs_init_pdfxtmd_cpdf( setLabel ,setMember ,beamID )
  character(*),intent(in) :: setLabel
  integer,intent(in),optional :: setMember,beamID
  integer :: idx,memberVal
  character(kind=c_char),allocatable :: cset(:)
  idx = 1
  if (present(beamID)) idx = max(1,min(size(pdfx_cpdf_state),beamID))
  memberVal = 0
  if (present(setMember)) memberVal = setMember
  call pdfx_ensure_cpdf_factory()
  if (.not. c_associated(pdfx_cpdf_factory)) then
    call pdfx_placeholder_notice('PDFxTMD cpdf factory missing',idx)
    return
  endif
  call pdfx_make_cstring(setLabel,cset)
  pdfx_cpdf_state(idx)%cpdfHandle = pdfx_mk_cpdf(pdfx_cpdf_factory,cset,memberVal)
  deallocate(cset)
  if (.not. c_associated(pdfx_cpdf_state(idx)%cpdfHandle)) then
    pdfx_cpdf_state(idx)%isReady = .false.
    call pdfx_placeholder_notice('mk_cpdf returned null',idx)
    return
  endif
  pdfx_cpdf_state(idx)%setLabel = trim(setLabel)
  pdfx_cpdf_state(idx)%member = memberVal
  pdfx_cpdf_state(idx)%isReady = .true.
  end subroutine

  function pdfvec_pdfxtmd( xx ,qq ,kTsq ,kf ,beamID ) result(list)
  intent(in) :: kf,xx,qq,kTsq
  integer,intent(in),optional :: beamID
  !(realknd2!) :: xx,qq,kTsq,list(-6:6)
  integer :: idx,kf
  list = 0
  idx = 1
  if (present(beamID)) idx = max(1,min(size(pdfx_state),beamID))
  list = pdfx_eval_tmd(idx,xx,kTsq,qq)
  end function

  function pdfvec_pdfxtmd_cpdf( xx ,qq ,beamID ) result(list)
  intent(in) :: xx,qq
  integer,intent(in),optional :: beamID
  !(realknd2!) :: xx,qq,list(-6:6)
  integer :: idx
  list = 0
  idx = 1
  if (present(beamID)) idx = max(1,min(size(pdfx_cpdf_state),beamID))
  list = pdfx_eval_cpdf(idx,xx,qq)
  end function

  function pdfunc_pdfxtmd( iParton ,xx ,qq ,kTsq ,kf ,beamID ) result(rslt)
  intent(in) :: iParton,kf,xx,qq,kTsq
  integer :: iParton,kf
  integer,intent(in),optional :: beamID
  !(realknd2!) :: xx,qq,kTsq,rslt,list(-6:6)
  list = pdfvec_pdfxtmd( xx ,qq * qq ,kTsq ,kf ,beamID )
  rslt = list(lgn2pdf(iParton))
  end function

  function pdfunc_pdfxtmd_cpdf( iParton ,xx ,qq ,beamID ) result(rslt)
  intent(in) :: iParton,xx,qq
  integer :: iParton
  integer,intent(in),optional :: beamID
  !(realknd2!) :: xx,qq,rslt,list(-6:6)
  list = pdfvec_pdfxtmd_cpdf( xx ,qq ,beamID )
  rslt = list(lgn2pdf(iParton))
  end function

  function pdfvec_pdfxtmdset( xx ,qq ,kTsq ,kf ,setLabel ,setMember ) result(list)
  intent(in) :: kf,setMember,xx,qq,kTsq
  character(*),intent(in) :: setLabel
  integer :: idx,kf,setMember
  !(realknd2!) :: xx,qq,kTsq,list(-6:6)
  idx = pdfx_match_state(setLabel,setMember,1)
  list = pdfx_eval_tmd(idx,xx,kTsq,qq * qq)
  end function

  subroutine pdfs_init_pdfxtmd_coupling( setLabel )
  character(*),intent(in) :: setLabel
  character(kind=c_char),allocatable :: cset(:)
  call pdfx_ensure_coupling_factory()
  if (.not. c_associated(pdfx_coupling_factory)) then
    call pdfx_placeholder_notice('PDFxTMD coupling factory missing',1)
    return
  endif
  call pdfx_make_cstring(setLabel,cset)
  pdfx_coupling_handle = pdfx_mk_coupling(pdfx_coupling_factory,cset)
  deallocate(cset)
  if (.not. c_associated(pdfx_coupling_handle)) then
    call pdfx_placeholder_notice('mk_coupling returned null',1)
    pdfx_use_coupling = .false.
    return
  endif
  pdfx_coupling_set = trim(setLabel)
  pdfx_use_coupling = .true.
  end subroutine

  function pdfunc_pdfxtmdset( iParton ,xx ,qq ,kTsq ,kf ,setLabel ,setMember ) result(rslt)
  intent(in) :: iParton,kf,setMember,xx,qq,kTsq
  integer :: iParton,kf,setMember
  character(*),intent(in) :: setLabel
  !(realknd2!) :: xx,qq,kTsq,rslt,list(-6:6)
  list = pdfvec_pdfxtmdset( xx ,qq ,kTsq ,kf ,setLabel ,setMember )
  rslt = list(lgn2pdf(iParton))
  end function

  subroutine pdfx_placeholder_notice( context ,beamID )
  character(*),intent(in) :: context
  integer,intent(in) :: beamID
  if (pdfx_warned) return
  if (errru.ge.0) then
    write(errru,*) 'INFO: PDFxTMD placeholder reached during ',trim(context)
    if (beamID>=1 .and. beamID<=size(pdfx_state)) then
      write(errru,*) '      Beam index: ',beamID,' set label: ',trim(pdfx_state(beamID)%setLabel)
    else
      write(errru,*) '      Beam index: ',beamID
    endif
    write(errru,*) '      Replace or fix PDFxTMD bindings.'
  endif
  pdfx_warned = .true.
  end subroutine

end module
