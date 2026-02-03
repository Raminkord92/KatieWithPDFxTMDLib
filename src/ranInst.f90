module katie_ranInst
  use avh_mathcnst
  use avh_iounits
  use avh_grid
  use avh_trivinp
  use avh_kaleu_mulcha, only: rng_interface
  use !(rngmodule!)
  implicit none
  private

  public :: noPdfs_init,noPdfs_gnrt,noPdfs_wght
  public :: noPdfs_collect,noPdfs_plot
  public :: ranInst_type


  type :: ranInst_type
    private
    type(grid_type) :: grdT,grdR,grdKTA,grdKTB
    !(realknd2!) :: Esq,Msq,a,b,c,normA,normB,xAmin,xAvol
    integer :: Nfinst,task,option,kinematics
    procedure(rng_interface),pointer,nopass :: rng
  contains
    procedure :: init_2
    procedure :: init_1
    generic :: init=>init_2,init_1
    procedure :: gnrt_AB
    procedure :: gnrt_A0
    procedure :: gnrt_0B
    procedure :: gnrt_00
    procedure :: gnrt_A
    procedure :: gnrt_0
    procedure :: adapt=>inst_adapt
    procedure :: plot=>inst_plot
    procedure :: dump=>inst_dump
  end type


  integer,save :: ostask
  type(grid_type),save :: phiA,phiB
  !(realknd2!),save :: Qz

contains
  

  subroutine init_2( obj ,rng ,task ,Nfinst ,EbeamA ,EbeamB &
                    ,Mass ,option ,xmin ,expon ,fileName ,readUnit )
  class(ranInst_type),intent(out) :: obj
  procedure(rng_interface) :: rng
  integer,intent(in) :: task(-2:-1)
  integer,intent(in) :: Nfinst,option,readUnit
  !(realknd2!),intent(in) :: EbeamA,EbeamB,Mass,xmin,expon
  character(*),intent(in) :: fileName
  optional :: Mass,option,xmin,expon,fileName,readUnit
  !(realknd2!) :: tmin
  integer :: iunit
!
  obj%rng => rng
  obj%Nfinst = Nfinst
  obj%Esq = 4*EbeamA*EbeamB
  obj%normA = log(1+EbeamA*EbeamA)
  obj%normB = log(1+EbeamB*EbeamB)
  obj%Msq = 125*125 ;if (present(Mass)) obj%Msq = Mass*Mass
  obj%option = 1  ;if (present(option)) obj%option = option 
  obj%task = 2*task(-2) + task(-1)
  obj%kinematics = 2
!
  if (present(fileName)) then
    iunit = deffile
    if (present(readUnit)) iunit = readUnit
    call obj%grdR%load( trim(fileName)//'r.instdump' ,iunit )
    call obj%grdT%load( trim(fileName)//'t.instdump' ,iunit )
    if (obj%task.eq.1.or.obj%task.eq.3) &
      call obj%grdKTA%load( trim(fileName)//'kTA.instdump' ,iunit )
    if (obj%task.eq.2.or.obj%task.eq.3) &
      call obj%grdKTB%load( trim(fileName)//'kTB.instdump' ,iunit )
  else
    call obj%grdR%init(nbatch=1000,nchmax=200,label='xAxBr')
    call obj%grdT%init(nbatch=1000,nchmax=200,label='xAxBt')
    call obj%grdKTA%init(nbatch=1000,nchmax=200,label='kTA')
    call obj%grdKTB%init(nbatch=1000,nchmax=200,label='kTB')
  endif
!
  select case (obj%option)
  case (1)
    tmin = 1d-16  ;if (present(xmin)) tmin = xmin**2
    obj%a = 1-log(tmin)
    obj%b = 1/obj%a
    obj%c = tmin*obj%a
  case (2)
    obj%a =-0.99d0  ;if (present(expon)) obj%a = expon
    if (obj%a.le.-1d0) obj%a =-0.99d0
    obj%c = 1/(1+obj%a)
  case default
    obj%a =-0.01d0  ;if (present(expon)) obj%a = expon
    if (obj%a.le.-1d0) obj%a =-0.99d0
    obj%c = 1/(1+obj%a)
  end select
  end subroutine
 

  subroutine init_1( obj ,rng ,task ,EbeamA ,EbeamB &
                    ,xmin ,xmax ,fileName ,readUnit )
  class(ranInst_type),intent(out) :: obj
  procedure(rng_interface) :: rng
  integer,intent(in) :: task(-2:-1)
  integer,intent(in) :: readUnit
  !(realknd2!),intent(in) :: EbeamA,EbeamB ,xmin,xmax
  character(*),intent(in) :: fileName
  optional :: xmin,xmax,fileName,readUnit
  !(realknd2!) :: tmin
  integer :: iunit
  obj%rng => rng
  obj%normA = log(1+EbeamA*EbeamA) 
  obj%normB = log(1+EbeamB*EbeamB) 
  obj%task = task(-1)
  obj%kinematics = 1
  obj%xAmin = 0           ;if (present(xmin)) obj%xAmin=xmin
  obj%xAvol = 1-obj%xAmin ;if (present(xmax)) obj%xAvol=xmax-obj%xAmin
  if (present(fileName)) then
    iunit = deffile
    if (present(readUnit)) iunit = readUnit
    call obj%grdR%load( trim(fileName)//'r.instdump' ,iunit )
    if (obj%task.eq.1) &
      call obj%grdKTA%load( trim(fileName)//'kTA.instdump' ,iunit )
  else
    call obj%grdR%init(nbatch=1000,nchmax=200,label='xA')
    call obj%grdKTA%init(nbatch=1000,nchmax=200,label='kTA')
  endif
  end subroutine
 

  subroutine inst_dump( obj ,fileName ,writeUnit )
  class(ranInst_type) :: obj
  character(*),intent(in) :: fileName
  integer,intent(in),optional :: writeUnit
  integer :: iunit
  iunit = deffile
  if (present(writeUnit)) iunit = writeUnit
  select case (obj%kinematics)
  case default
    call obj%grdR%dump( iunit ,trim(fileName)//'r.instdump' )
    call obj%grdT%dump( iunit ,trim(fileName)//'t.instdump' )
    if (obj%task.eq.1.or.obj%task.eq.3) &
      call obj%grdKTA%dump( iunit ,trim(fileName)//'kTA.instdump' )
    if (obj%task.eq.2.or.obj%task.eq.3) &
      call obj%grdKTB%dump( iunit ,trim(fileName)//'kTB.instdump' )
  case (1)
    call obj%grdR%dump( iunit ,trim(fileName)//'r.instdump' )
    if (obj%task.eq.1) &
      call obj%grdKTA%dump( iunit ,trim(fileName)//'kTA.instdump' )
  end select
  end subroutine


  subroutine xAxB_gnrt( obj ,ww ,xA,xB )
  class(ranInst_type) :: obj
  !(realknd2!),intent(  out) :: ww,xA,xB
  !(realknd2!) :: rr,tt,hh
  tt = obj%grdT%gnrt( obj%rng() )
  rr = obj%grdR%gnrt( obj%rng() )
  ww = obj%grdT%wght() * obj%grdR%wght()
  select case (obj%option)
  case (1) 
! Generate  rr=ln(xA)/(ln(xA)+ln(xB))  and  tt=xA*xB
! tt in [0,1] following  theta(tt<tmin)/tmin + theta(tt>tmin)/tt
    if (tt.le.obj%b) then
      tt = tt*obj%c
      ww = ww*obj%c/tt
      hh = log(tt)
    else
      ww = ww*obj%a
      hh = (tt-1)*obj%a
    endif
    xA = exp( hh*(1-rr) )
    xB = exp( hh*(  rr) )
    ww =-ww*hh*xA*xB
  case (2)
! Generate  rr=ln(xA)/(ln(xA)+ln(xB))  and  tt=xA*xB
! tt in [0,1] following x^(obj%a)
    hh = tt**(-obj%a*obj%c)
    tt = tt*hh
    ww = ww*hh*obj%c*(-log(tt))
    xB = tt**rr
    xA = tt/xB
  case default
! xA and xB independent following x^(obj%a)
    hh = tt**(-obj%a*obj%c)
      ww = ww*hh*obj%c
      xA = tt*hh
    hh = rr**(-obj%a*obj%c)
      ww = ww*hh*obj%c
      xB = rr*hh
  end select
  end subroutine
 

  subroutine gnrt_AB( obj ,weight ,xA,kTA ,xB,kTB )
  class(ranInst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xA,kTA(2) ,xB,kTB(2)
  !(realknd2!) :: rhoA,rhoB,wght,ww,QA,QB,phi,cp,sp,ct,st,rA(1:2),rB(1:2)
  weight = 0
  ww = 1
  call xAxB_gnrt( obj ,wght ,rhoA,rhoB ) ;ww=ww*wght
  xA = rhoA
  xB = rhoB
  QA = exp( obj%normA * obj%grdKTA%gnrt(obj%rng()) )
  ww = ww * obj%normA*QA * obj%grdKTA%wght()
  QA = sqrt(QA-1)
  QB = exp( obj%normB * obj%grdKTB%gnrt(obj%rng()) )
  ww = ww * obj%normB*QB * obj%grdKTB%wght()
  QB = sqrt(QB-1)
  if (obj%Nfinst.gt.1) then
    phi = r2PI*obj%rng()
    kTA = [QA*cos(phi),QA*sin(phi)]
    ww = ww * r1PI
    phi = r2PI*obj%rng()
    kTB = [QB*cos(phi),QB*sin(phi)]
    ww = ww * r1PI
  else
    ct = xA*xB*obj%Esq-obj%Msq
    if (ct.le.rZRO) return
    ct = (ct - QA*QA - QB*QB)/(2*QA*QB)
    if (abs(ct).gt.rONE) return
    st = sqrt(1-ct*ct)
    rA(1:2) = [ rZRO,   QA]
    rB(1:2) = [st*QB,ct*QB]
    phi = r2PI*obj%rng()
    cp = cos(phi)
    sp = sin(phi)
    kTA(1:2) = [         sp*rA(2),          cp*rA(2)]
    kTB(1:2) = [cp*rB(1)+sp*rB(2),-sp*rB(1)+cp*rB(2)]
    ww = ww * r1PI/(st*QA*QB) / 2
  endif
  weight = ww
  end subroutine


  subroutine gnrt_A0( obj ,weight ,xA,kTA ,xB,kTB )
  class(ranInst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xA,kTA(2) ,xB,kTB(2)
  !(realknd2!) :: rhoA,rhoB,wght,ww,QA,phi
  weight = 0
  ww = 1
  call xAxB_gnrt( obj ,wght ,rhoA,rhoB ) ;ww=ww*wght
  xA = rhoA
  xB = rhoB
  if (obj%Nfinst.gt.1) then
    QA = exp( obj%normA * obj%grdKTA%gnrt(obj%rng()) )
    ww = ww * r1PI * obj%normA*QA * obj%grdKTA%wght()
    QA = sqrt(QA-1)
  else
    QA = xA*xB*obj%Esq-obj%Msq
    if (QA.le.rZRO) return
    QA = sqrt(QA)
    ww = ww * r1PI!*QA
  endif
  phi = r2PI*obj%rng()
  kTA = [QA*cos(phi),QA*sin(phi)]
  kTB = 0
  weight = ww
  end subroutine


  subroutine gnrt_0B( obj ,weight ,xA,kTA ,xB,kTB )
  class(ranInst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xA,kTA(2) ,xB,kTB(2)
  !(realknd2!) :: rhoA,rhoB,wght,ww,QB,phi
  weight = 0
  ww = 1
  call xAxB_gnrt( obj ,wght ,rhoA,rhoB ) ;ww=ww*wght
  xA = rhoA
  xB = rhoB
  if (obj%Nfinst.gt.1) then
    QB = exp( obj%normB * obj%grdKTB%gnrt(obj%rng()) )
    ww = ww * r1PI * obj%normB*QB * obj%grdKTB%wght()
    QB = sqrt(QB-1)
  else
    QB = xA*xB*obj%Esq-obj%Msq
    if (QB.le.rZRO) return
    QB = sqrt(QB)
    ww = ww * r1PI!*QB
  endif
  phi = r2PI*obj%rng()
  kTB = [QB*cos(phi),QB*sin(phi)]
  kTA = 0
  weight = ww
  end subroutine


  subroutine gnrt_00( obj ,weight ,xA,kTA ,xB,kTB )
  class(ranInst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xA,kTA(2) ,xB,kTB(2)
  !(realknd2!) :: rhoA,rhoB,wght,ww,phi,tt
  weight = 0
  ww = 1
  if (obj%Nfinst.gt.1) then
    call xAxB_gnrt( obj ,wght ,rhoA,rhoB ) ;ww=ww*wght
    xA = rhoA
    xB = rhoB
  else
    tt = obj%Msq/obj%Esq
    xB = tt**obj%grdT%gnrt( obj%rng() )
    xA = tt/xB
    ww =-log(tt)*obj%grdT%wght()/obj%Esq
  endif
  kTA = 0
  kTB = 0
  weight = ww
  end subroutine


  subroutine gnrt_A( obj ,weight ,xA,kTA )
  class(ranInst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xA,kTA(2)
  !(realknd2!) :: QA,phi,wght
  xA = obj%grdR%gnrt( obj%rng() )
  weight = obj%grdR%wght()
  xA = obj%xAvol*xA + obj%xAmin
  weight = weight*obj%xAvol
  QA = exp( obj%normA * obj%grdKTA%gnrt(obj%rng()) )
  weight = weight * r1PI * obj%normA*QA * obj%grdKTA%wght()
  QA = sqrt(QA-1)
  phi = r2PI*obj%rng()
  kTA = [QA*cos(phi),QA*sin(phi)]
  end subroutine

  subroutine gnrt_0( obj ,weight ,xA,kTA )
  class(ranInst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xA,kTA(2)
  xA = obj%grdR%gnrt( obj%rng() )
  weight = obj%grdR%wght()
  xA = obj%xAvol*xA + obj%xAmin
  weight = weight*obj%xAvol
  kTA = 0
  end subroutine


  subroutine inst_adapt( obj ,weight )
  class(ranInst_type) :: obj
  !(realknd2!),intent(in) :: weight
  if (weight.le.rZRO) return
  select case (obj%kinematics)
  case default
    select case (obj%task)
    case (3)
      call obj%grdT%collect( weight )
      call obj%grdR%collect( weight )
      call obj%grdKTA%collect( weight )
      call obj%grdKTB%collect( weight )
    case (2)
      call obj%grdT%collect( weight )
      call obj%grdR%collect( weight )
      if (obj%Nfinst.gt.1) call obj%grdKTB%collect( weight )
    case (1)
      call obj%grdT%collect( weight )
      call obj%grdR%collect( weight )
      if (obj%Nfinst.gt.1) call obj%grdKTA%collect( weight )
    case (0)
      call obj%grdT%collect( weight )
      if (obj%Nfinst.gt.1) call obj%grdR%collect( weight )
    end select
  case (1)
    select case (obj%task)
    case (1)
      call obj%grdR%collect( weight )
      call obj%grdKTA%collect( weight )
    case (0)
      call obj%grdR%collect( weight )
    end select
  end select
  end subroutine


  subroutine inst_plot( obj ,iunit )
  class(ranInst_type) :: obj
  integer,intent(in) ::iunit
  select case (obj%kinematics)
  case default
    select case (obj%task)
    case (3)
      call obj%grdT%plot( iunit )
      call obj%grdR%plot( iunit )
      call obj%grdKTA%plot( iunit )
      call obj%grdKTB%plot( iunit )
    case (2)
      call obj%grdT%plot( iunit )
      call obj%grdR%plot( iunit )
      if (obj%Nfinst.gt.1) call obj%grdKTB%plot( iunit )
    case (1)
      call obj%grdT%plot( iunit )
      call obj%grdR%plot( iunit )
      if (obj%Nfinst.gt.1) call obj%grdKTA%plot( iunit )
    case (0)
      call obj%grdT%plot( iunit )
      if (obj%Nfinst.gt.1) call obj%grdR%plot( iunit )
    end select
  case (1)
    select case (obj%task)
    case (1)
      call obj%grdR%plot( iunit )
      call obj%grdKTA%plot( iunit )
    case (0)
      call obj%grdR%plot( iunit )
    end select
  end select
  end subroutine

 
  subroutine noPdfs_init(ctask,Q0)
  integer,intent(in) :: ctask(-2:-1)
  !(realknd2!),intent(in),optional :: Q0
  Qz = 1.0d0 ;if (present(Q0)) Qz=Q0
  ostask = 2*ctask(-2)+ctask(-1)
  select case (ostask)
  case (3)
    call phiA%init(nbatch=1000,nchmax=200,label='kAphi')
    call phiB%init(nbatch=1000,nchmax=200,label='kBphi')
  case (2)
    call phiB%init(nbatch=1000,nchmax=200,label='kBphi')
  case (1)
    call phiA%init(nbatch=1000,nchmax=200,label='kAphi')
  case (0)
  case default
    if (errru.ge.0) write(errru,*) 'ERROR in noPdfs: task not defined'
    stop
  end select
  end subroutine

  subroutine noPdfs_gnrt(kTA,kTB)
  !(realknd2!),intent(out) :: kTA(1:2),kTB(1:2)
  !(realknd2!) :: phi
  select case (ostask)
  case (3)
    !(rngenerator phi !)
    phi=phiA%gnrt(phi)*r2PI ;kTA=[Qz*cos(phi),Qz*sin(phi)]
    !(rngenerator phi !)
    phi=phiB%gnrt(phi)*r2PI ;kTB=[Qz*cos(phi),Qz*sin(phi)]
  case (2)
    kTA = 0
    !(rngenerator phi !)
    phi=phiB%gnrt(phi)*r2PI ;kTB=[Qz*cos(phi),Qz*sin(phi)]
  case (1)
    !(rngenerator phi !)
    phi=phiA%gnrt(phi)*r2PI ;kTA=[Qz*cos(phi),Qz*sin(phi)]
    kTB = 0
  case (0)
    kTA = 0
    kTB = 0
  end select
  end subroutine

  function noPdfs_wght() result(rslt)
  !(realknd2!) :: rslt
  select case (ostask)
  case (3) ;rslt = phiA%wght()*phiB%wght()
  case (2) ;rslt = phiB%wght()
  case (1) ;rslt = phiA%wght()
  case (0) ;rslt = 1
  end select
  end function

  subroutine noPdfs_collect(weight)
  !(realknd2!),intent(in) :: weight
  select case (ostask)
  case (3) 
    call phiA%collect( weight )
    call phiB%collect( weight )
  case (2) 
    call phiB%collect( weight )
  case (1)
    call phiA%collect( weight )
  end select
  end subroutine

  subroutine noPdfs_plot( iunit )
  integer,intent(in) ::iunit
  select case (ostask)
  case (3)
    call phiA%plot( iunit )
    call phiB%plot( iunit )
  case (2)
    call phiB%plot( iunit )
  case (1)
    call phiA%plot( iunit )
  end select
  end subroutine


end module


