module katie_DISinst
  use avh_mathcnst
  use avh_iounits
  use avh_grid
  use avh_kaleu_mulcha, only: rng_interface
  use !(rngmodule!)
  implicit none
  private

  public :: DISinst_type

  type :: DISinst_type
    private
    type(grid_type) :: grid_kT,grid_Msq,grid_xB,grid_Qsq
    !(realknd2!) :: Eh,Ee,xBmin,xBvol,QsqMin,QsqVol
    integer :: Nfinst,task
    procedure(rng_interface),pointer,nopass :: rng
  contains
    procedure :: init
    procedure :: dump
    procedure :: gnrt
    procedure :: adapt
    procedure :: plot
  end type


contains
  
  subroutine init( obj ,rng ,task ,Ehadr,Eelec ,Nfinst &
                  ,xBmin,xBmax ,QsqMin,QsqMax ,pathFile ,readUnit )
  class(DISinst_type) :: obj
  !(realknd2!),intent(in) :: Ehadr,Eelec ,xBmin,xBmax ,QsqMin,QsqMax
  procedure(rng_interface) :: rng
  integer,intent(in) :: Nfinst,task,readUnit
  character(*),intent(in) :: pathFile
  optional :: pathFile,readUnit
  integer :: iunit
  if (rZRO.gt.xBmin.or.xBmax.gt.rONE) then
    write(*,*) 'ERROR: incompatible values requested for xB'
    stop
  endif
  if (rZRO.gt.QsqMin.or.QsqMax.gt.4*Ehadr*Eelec) then
    write(*,*) 'ERROR: incompatible values requested for Qsq'
    stop
  endif
  obj%Eh = Ehadr
  obj%Ee = Eelec
  obj%rng => rng
  obj%xBmin = xBmin
  obj%xBvol = xBmax-xBmin
  obj%QsqMin = QsqMin
  obj%QsqVol = QsqMax-QsqMin
  obj%Nfinst = Nfinst
  obj%task = task
  if (present(pathFile)) then
    iunit = deffile
    if (present(readUnit)) iunit = readUnit
    if (obj%task.ne.0)   call obj%grid_kT%load( trim(pathFile)//'kT.instdump' ,iunit )
    if (obj%Nfinst.gt.1) call obj%grid_Msq%load( trim(pathFile)//'Msq.instdump' ,iunit )
    if (obj%xBvol.gt.rZRO)  call obj%grid_xB%load(  trim(pathFile)//'xB.instdump'  ,iunit )
    if (obj%QsqVol.gt.rZRO) call obj%grid_Qsq%load( trim(pathFile)//'Qsq.instdump' ,iunit )
  else
    if (obj%task.ne.0)   call obj%grid_kT%init( nbatch=1000 ,nchmax=200 ,label='kT' )
    if (obj%Nfinst.gt.1) call obj%grid_Msq%init( nbatch=1000 ,nchmax=200 ,label='Msq' )
    if (obj%xBvol.gt.rZRO)  call obj%grid_xB%init(  nbatch=1000 ,nchmax=200 ,label='xB'  )
    if (obj%QsqVol.gt.rZRO) call obj%grid_Qsq%init( nbatch=1000 ,nchmax=200 ,label='Qsq' )
  endif
  end subroutine

  subroutine dump( obj ,pathFile ,writeUnit )
  class(DISinst_type) :: obj
  character(*),intent(in) :: pathFile
  integer,intent(in),optional :: writeUnit
  integer :: iunit
  iunit = deffile
  if (present(writeUnit)) iunit = writeUnit
  if (obj%task.ne.0)   call obj%grid_kT%dump( iunit ,trim(pathFile)//'kT.instdump' )
  if (obj%Nfinst.gt.1) call obj%grid_Msq%dump( iunit ,trim(pathFile)//'Msq.instdump' )
  if (obj%xBvol.gt.rZRO)  call obj%grid_xB%dump(  iunit ,trim(pathFile)//'xB.instdump'  )
  if (obj%QsqVol.gt.rZRO) call obj%grid_Qsq%dump( iunit ,trim(pathFile)//'Qsq.instdump' )
  end subroutine

  subroutine plot( obj ,pathFile ,writeUnit )
  class(DISinst_type) :: obj
  character(*),intent(in) :: pathFile
  integer,intent(in),optional :: writeUnit
  integer :: iunit
  iunit = deffile
  if (present(writeUnit)) iunit = writeUnit
  if (obj%task.ne.0)   call obj%grid_kT%plot( iunit ,trim(pathFile)//'kT.plot' )
  if (obj%Nfinst.gt.1) call obj%grid_Msq%plot( iunit ,trim(pathFile)//'Msq.plot' )
  if (obj%xBvol.gt.rZRO)  call obj%grid_xB%plot(  iunit ,trim(pathFile)//'xB.plot'  )
  if (obj%QsqVol.gt.rZRO) call obj%grid_Qsq%plot( iunit ,trim(pathFile)//'Qsq.plot' )
  end subroutine

  subroutine adapt( obj ,weight )
  class(DISinst_type) :: obj
  !(realknd2!),intent(in) :: weight
  if (weight.le.rZRO) return
  if (obj%task.ne.0)   call obj%grid_kT%collect( weight )
  if (obj%Nfinst.gt.1) call obj%grid_Msq%collect( weight )
  if (obj%xBvol.gt.rZRO)  call obj%grid_xB%collect( weight )
  if (obj%QsqVol.gt.rZRO) call obj%grid_Qsq%collect( weight )
  end subroutine

  subroutine gnrt( obj ,weight ,xB,Qsq,yy ,xx,kTsq ,qMOM,kMOM )
  class(DISinst_type) :: obj
  !(realknd2!),intent(out) :: weight ,xB,Qsq,yy ,xx,kTsq ,kMOM(0:3),qMOM(0:3)
  !(realknd2!) :: vv,Msq,kTsqMx,phiKt,phiQt,aa,bb,Qabs,kTabs,kT(3),qT(2),MsqMax
  !(realknd2!) :: wght
  weight = 0
  wght = 1
!
  xB = obj%xBmin
  if (obj%xBvol.gt.rZRO) then
    xB = xB + obj%xBvol*obj%grid_xB%gnrt(obj%rng())
    wght = wght * obj%xBvol * obj%grid_xB%wght()
  endif
  Qsq = obj%QsqMin
  if (obj%QsqVol.gt.rZRO) then
    Qsq = Qsq + obj%QsqVol*obj%grid_Qsq%gnrt(obj%rng())
    wght = wght * obj%QsqVol * obj%grid_Qsq%wght()
  endif
  Qabs = sqrt(Qsq)
  yy = Qsq/(xB*4*obj%Eh*obj%Ee)
  vv = 1-yy
  if (vv.lt.rZRO) return
  vv = sqrt(vv)
  Msq = 0
  if (obj%Nfinst.gt.1) then
    MsqMax = (1/xB-yy)*Qsq
    Msq = Msq + MsqMax*obj%grid_Msq%gnrt(obj%rng())
    wght = wght * MsqMax * obj%grid_Msq%wght()
  endif
!
  if (obj%task.eq.0) then
    kTsq = rZRO
    kTabs = 0
    phiKt = 0
    kT(1:2) = 0
    wght = wght*2 ! to correct the 1/32 to 1/16 below for on-shell
  else
    kTsqMx = 1/xB-yy-Msq/Qsq
    if (kTsqMx.lt.rZRO) return
    kTsqMx = Qsq*( vv + sqrt(kTsqMx) )**2
    !OLD kTsq = kTsqMx*obj%grid_kT%gnrt(obj%rng())
    !OLD wght = wght * kTsqMx * obj%grid_kT%wght()
    kTsqMx = log(1+kTsqMx)                               !NEW
    kTsq = exp( kTsqMx*obj%grid_kT%gnrt(obj%rng()) )-1   !NEW
    wght = wght * kTsqMx * (1+kTsq) * obj%grid_kT%wght() !NEW
    phiKt = r2PI*obj%rng()
    wght = wght * r2PI
    kTabs = sqrt(kTsq)
    kT(1:2) = kTabs*[sin(phiKt),cos(phiKt)]
  endif
!
  phiQt = r2PI*obj%rng()
  wght = wght * r2PI
  xx = xB*( 1+(kTsq+Msq)/Qsq + 2*cos(phiQt)*vv*kTabs/Qabs )
  qT(1:2) = vv*Qabs*[sin(phiKt+phiQt),cos(phiKt+phiQt)]
!
  aa = Qsq/(4*xB*obj%Eh)
  bb =-Qsq/(4*obj%Ee)
  qMOM = [   aa+bb ,qT(1),qT(2),   bb-aa ]
  kMOM = [xx*obj%Eh,kT(1),kT(2),xx*obj%Eh]
! 
  wght = wght/(32*xB*obj%Eh*obj%Ee)
  weight = wght
  end subroutine


end module


