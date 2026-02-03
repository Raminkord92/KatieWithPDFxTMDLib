module katie_histogramtools
  implicit none
  private
  public :: histo_1d_type ,histo_2d_type
  public :: mass,pTrans,ETrans,rapidity,pseudoRap,phi,deltaPhi,deltaR
  public :: theta,angle,breit0_type,breit_type
  public :: unity,numeral,r1PI,r2PI,sort_big2small,sort_small2big

  type :: histo_1d_type
    real(kind(1d0)),allocatable :: left(:),rght(:),sumw(:,:),Nev
    integer :: Nbins
  contains
    procedure :: init_1d_a
    procedure :: init_1d_b
    procedure :: init_1d_c
    procedure :: init_1d_d
    procedure :: init_1d_e
    generic :: init => init_1d_a,init_1d_b,init_1d_c,init_1d_d,init_1d_e
    procedure :: collect => collect_1d
    procedure :: write => write_1d
  end type

  type :: histo_2d_type
    real(kind(1d0)),allocatable :: xLow(:),xUpp(:),yLow(:),yUpp(:),sumw(:,:,:),Nev
    integer :: NxBins,NyBins
  contains
    procedure :: init_2d_a
    procedure :: init_2d_b
    procedure :: init_2d_c
    procedure :: init_2d_d
    procedure :: init_2d_e
    generic :: init => init_2d_a,init_2d_b,init_2d_c,init_2d_d,init_2d_e
    procedure :: collect => collect_2d
    procedure :: write => write_2d
  end type

  type :: breit0_type
    real(kind(1d0)) :: rr(0:3),mm,stsp,stcp,ct,st,sp,cp,ctsp,ctcp,sgn
  contains
    procedure :: init=>breit0_init
    procedure :: act=>breit0_act
    procedure :: act_ev=>breit0_act_ev
  end type

  type :: breit_type
    real(kind(1d0)) :: cp,sp,c00,c01,c03,c10,c13,c30,c31,c33
  contains
    procedure :: init=>breit_init
    procedure :: act=>breit_act
    procedure :: act_ev=>breit_act_ev
  end type

  integer,parameter :: unity(0:9) = [0,1,2,3,4,5,6,7,8,9]
  character(1),parameter :: numeral(0:9) = ['0','1','2','3','4','5','6','7','8','9']
  real(kind(1d0)),parameter :: r1PI = 3.1415926535897932384626433832795d0
  real(kind(1d0)),parameter :: r2PI = 2*r1PI

contains


  subroutine init_1d_a( obj ,left ,right ,Nbins )
  class(histo_1d_type) :: obj
  real(kind(1d0)),intent(in) :: left,right
  integer,intent(in) :: Nbins
  integer :: ii
  real(kind(1d0)) :: step
  if (allocated(obj%left)) deallocate(obj%left)
  if (allocated(obj%rght)) deallocate(obj%rght)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  allocate( obj%left(1:Nbins) ,obj%rght(1:Nbins) ,obj%sumw(1:2,1:Nbins) )
  obj%Nbins = Nbins
  step = (right-left)/Nbins
  do ii=1,Nbins
    obj%left(ii) = left + (ii-1)*step
    obj%rght(ii) = left + ii*step
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_1d_e( obj ,left ,right ,eps ,Nbins )
! The right end of the last bin is          left+(right-left)*eps^0 ,
! the right end of the next-to-last bin is  left+(right-left)*eps^1 ,
! the right end of the NN-last bin is       left+(right-left)*eps^2 ,
! etc.
  class(histo_1d_type) :: obj
  real(kind(1d0)),intent(in) :: left,right,eps
  integer,intent(in) :: Nbins
  integer :: ii
  real(kind(1d0)) :: step
  if (allocated(obj%left)) deallocate(obj%left)
  if (allocated(obj%rght)) deallocate(obj%rght)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  allocate( obj%left(1:Nbins) ,obj%rght(1:Nbins) ,obj%sumw(1:2,1:Nbins) )
  obj%Nbins = Nbins
  step = (right-left)
  do ii=0,Nbins-1
    obj%rght(Nbins-ii) = left + step
    step = step*eps
    obj%left(Nbins-ii) = left + step
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_1d_b( obj ,bins )
  class(histo_1d_type) :: obj
  real(kind(1d0)),intent(in) :: bins(:)
  integer :: ii
  if (allocated(obj%left)) deallocate(obj%left)
  if (allocated(obj%rght)) deallocate(obj%rght)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  obj%Nbins = size(bins)/2
  allocate(obj%left(1:obj%Nbins),obj%rght(1:obj%Nbins),obj%sumw(1:2,1:obj%Nbins))
  do ii=1,obj%Nbins
    obj%left(ii) = bins(2*ii-1)
    obj%rght(ii) = bins(2*ii  )
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_1d_c( obj ,bins )
  class(histo_1d_type) :: obj
  real(kind(1e0)),intent(in) :: bins(:)
  integer :: ii
  if (allocated(obj%left)) deallocate(obj%left)
  if (allocated(obj%rght)) deallocate(obj%rght)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  obj%Nbins = size(bins)/2
  allocate(obj%left(1:obj%Nbins),obj%rght(1:obj%Nbins),obj%sumw(1:2,1:obj%Nbins))
  do ii=1,obj%Nbins
    obj%left(ii) = bins(2*ii-1)
    obj%rght(ii) = bins(2*ii  )
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_1d_d( obj ,bins )
  class(histo_1d_type) :: obj
  integer,intent(in) :: bins(:)
  integer :: ii
  if (allocated(obj%left)) deallocate(obj%left)
  if (allocated(obj%rght)) deallocate(obj%rght)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  obj%Nbins = size(bins)/2
  allocate(obj%left(1:obj%Nbins),obj%rght(1:obj%Nbins),obj%sumw(1:2,1:obj%Nbins))
  do ii=1,obj%Nbins
    obj%left(ii) = bins(2*ii-1)
    obj%rght(ii) = bins(2*ii  )
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine


  subroutine collect_1d( obj ,val ,weight )
  class(histo_1d_type) :: obj
  real(kind(1d0)),intent(in) :: val,weight
  integer :: ii,i0,i1
  obj%Nev = obj%Nev+1
  if (weight.eq.0d0) return
  if (val.lt.obj%left(1)) return
  if (val.ge.obj%rght(obj%Nbins)) return
  i0 = 1
  i1 = obj%Nbins
  do while (i0.ne.i1)
    ii = (i0+i1)/2
    if (val.lt.obj%left(ii)) then
      i1 = ii-1
    elseif (val.lt.obj%rght(ii)) then
      i1 = ii
      exit
    elseif (val.ge.obj%rght(ii)) then
      i0 = ii+1
    endif
  enddo
  if (val.lt.obj%left(i1).or.obj%rght(i1).le.val) return
  obj%sumw(1:2,i1) = obj%sumw(1:2,i1) + [weight,weight*weight]
  end subroutine


  subroutine write_1d( obj ,filename ,frmt ,factor ,writeUnit )
  class(histo_1d_type) :: obj
  character(*),intent(in) :: filename
  character(*),intent(in) :: frmt
  real(kind(1d0)),intent(in) :: factor
  integer,intent(in) :: writeUnit
  optional :: frmt,factor,writeUnit
  integer :: ii,wUnit
  real(kind(1d0)) :: h0,h1,h2
  h0 = obj%Nev
  wUnit=59;if(present(writeUnit))wUnit=writeUnit
  open(wUnit,file=trim(filename),status='new')
  !This is f2008: open(newunit=wUnit,file=trim(filename),status='new')
  select case (frmt)
  case default
    do ii=1,obj%Nbins
      call compute
      write(wUnit,'(4e16.8)') obj%left(ii),obj%rght(ii),h1,h2
    enddo
  case ('noBinNorm')
    do ii=1,obj%Nbins
      call compute
      h1 = h1*(obj%rght(ii)-obj%left(ii))
      h2 = h2*(obj%rght(ii)-obj%left(ii))
      write(wUnit,'(3e16.8)') obj%left(ii),h1,h2
      write(wUnit,'(3e16.8)') obj%rght(ii),h1,h2
    enddo
  case ('boxes')
    do ii=1,obj%Nbins
      call compute
      write(wUnit,'(3e16.8)') obj%left(ii),0d0,0d0
      write(wUnit,'(3e16.8)') obj%left(ii),h1,h2
      write(wUnit,'(3e16.8)') obj%rght(ii),h1,h2
      write(wUnit,'(3e16.8)') obj%rght(ii),0d0,0d0
    enddo
  case ('plateaux')
    do ii=1,obj%Nbins
      call compute
      write(wUnit,'(3e16.8)') obj%left(ii),h1,h2
      write(wUnit,'(3e16.8)') obj%rght(ii),h1,h2
    enddo
  end select
  close(wUnit)
  contains
    subroutine compute
    if (h0.le.0) then
      h1 = 0
      h2 = 0
    else
      h1 = obj%sumw(1,ii)/h0
      if (h0.le.1) then
        h2 = abs(h1)
      else
        h2 = (obj%sumw(2,ii)/h0 - h1*h1)/(h0-1)
        if (h2.lt.0) then
          h2 = abs(h1)
        else
          h2 = sqrt(h2)
        endif
      endif
    endif
    h1 = h1/(obj%rght(ii)-obj%left(ii))
    h2 = h2/(obj%rght(ii)-obj%left(ii))
    if (present(factor)) then
      h1 = h1*factor
      h2 = h2*factor
    endif
    end subroutine
  end subroutine


  subroutine init_2d_a( obj ,xLow,xUpp,NxBins ,yLow,yUpp,NyBins )
  class(histo_2d_type) :: obj
  real(kind(1d0)),intent(in) :: xLow,xUpp,yLow,yUpp
  integer,intent(in) :: NxBins,NyBins
  integer :: ii
  real(kind(1d0)) :: step
  if (allocated(obj%xLow)) deallocate(obj%xLow)
  if (allocated(obj%xUpp)) deallocate(obj%xUpp)
  if (allocated(obj%yLow)) deallocate(obj%yLow)
  if (allocated(obj%yUpp)) deallocate(obj%yUpp)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  allocate( obj%xLow(1:NxBins) ,obj%xUpp(1:NxBins) )
  allocate( obj%yLow(1:NyBins) ,obj%yUpp(1:NyBins) )
  allocate( obj%sumw(1:2,1:NxBins,1:NyBins) )
  obj%NxBins = NxBins
  obj%NyBins = NyBins
  step = (xUpp-xLow)/NxBins
  do ii=1,NxBins
    obj%xLow(ii) = xLow + (ii-1)*step
    obj%xUpp(ii) = xLow + ii*step
  enddo
  step = (yUpp-yLow)/NyBins
  do ii=1,NyBins
    obj%yLow(ii) = yLow + (ii-1)*step
    obj%yUpp(ii) = yLow + ii*step
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_2d_e( obj ,xLow,xUpp,xEps,NxBins ,yLow,yUpp,yEps,NyBins )
  class(histo_2d_type) :: obj
  real(kind(1d0)),intent(in) :: xLow,xUpp,xEps,yLow,yUpp,yEps
  integer,intent(in) :: NxBins,NyBins
  integer :: ii
  real(kind(1d0)) :: step
  if (allocated(obj%xLow)) deallocate(obj%xLow)
  if (allocated(obj%xUpp)) deallocate(obj%xUpp)
  if (allocated(obj%yLow)) deallocate(obj%yLow)
  if (allocated(obj%yUpp)) deallocate(obj%yUpp)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  allocate( obj%xLow(1:NxBins) ,obj%xUpp(1:NxBins) )
  allocate( obj%yLow(1:NyBins) ,obj%yUpp(1:NyBins) )
  allocate( obj%sumw(1:2,1:NxBins,1:NyBins) )
  obj%NxBins = NxBins
  obj%NyBins = NyBins
  step = xUpp-xLow
  do ii=0,NxBins-1
    obj%xUpp(NxBins-ii) = xLow + step
    step = step*xEps
    obj%xLow(NxBins-ii) = xLow + step
  enddo
  step = yUpp-yLow
  do ii=0,NyBins-1
    obj%yUpp(NyBins-ii) = yLow + step
    step = step*yEps
    obj%yLow(NyBins-ii) = yLow + step
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_2d_b( obj ,xBins ,yBins )
  class(histo_2d_type) :: obj
  real(kind(1d0)),intent(in) :: xBins(:),yBins(:)
  integer :: ii
  if (allocated(obj%xLow)) deallocate(obj%xLow)
  if (allocated(obj%xUpp)) deallocate(obj%xUpp)
  if (allocated(obj%yLow)) deallocate(obj%yLow)
  if (allocated(obj%yUpp)) deallocate(obj%yUpp)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  obj%NxBins = size(xBins)/2
  obj%NyBins = size(yBins)/2
  allocate( obj%xLow(1:obj%NxBins) ,obj%xUpp(1:obj%NxBins) )
  allocate( obj%yLow(1:obj%NyBins) ,obj%yUpp(1:obj%NyBins) )
  allocate( obj%sumw(1:2,1:obj%NxBins,1:obj%NyBins) )
  do ii=1,obj%NxBins
    obj%xLow(ii) = xBins(2*ii-1)
    obj%xUpp(ii) = xBins(2*ii  )
  enddo
  do ii=1,obj%NyBins
    obj%yLow(ii) = yBins(2*ii-1)
    obj%yUpp(ii) = yBins(2*ii  )
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_2d_c( obj ,xBins ,yBins )
  class(histo_2d_type) :: obj
  real(kind(1e0)),intent(in) :: xBins(:),yBins(:)
  integer :: ii
  if (allocated(obj%xLow)) deallocate(obj%xLow)
  if (allocated(obj%xUpp)) deallocate(obj%xUpp)
  if (allocated(obj%yLow)) deallocate(obj%yLow)
  if (allocated(obj%yUpp)) deallocate(obj%yUpp)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  obj%NxBins = size(xBins)/2
  obj%NyBins = size(yBins)/2
  allocate( obj%xLow(1:obj%NxBins) ,obj%xUpp(1:obj%NxBins) )
  allocate( obj%yLow(1:obj%NyBins) ,obj%yUpp(1:obj%NyBins) )
  allocate( obj%sumw(1:2,1:obj%NxBins,1:obj%NyBins) )
  do ii=1,obj%NxBins
    obj%xLow(ii) = xBins(2*ii-1)
    obj%xUpp(ii) = xBins(2*ii  )
  enddo
  do ii=1,obj%NyBins
    obj%yLow(ii) = yBins(2*ii-1)
    obj%yUpp(ii) = yBins(2*ii  )
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine

  subroutine init_2d_d( obj ,xBins ,yBins )
  class(histo_2d_type) :: obj
  integer,intent(in) :: xBins(:),yBins(:)
  integer :: ii
  if (allocated(obj%xLow)) deallocate(obj%xLow)
  if (allocated(obj%xUpp)) deallocate(obj%xUpp)
  if (allocated(obj%yLow)) deallocate(obj%yLow)
  if (allocated(obj%yUpp)) deallocate(obj%yUpp)
  if (allocated(obj%sumw)) deallocate(obj%sumw)
  obj%NxBins = size(xBins)/2
  obj%NyBins = size(yBins)/2
  allocate( obj%xLow(1:obj%NxBins) ,obj%xUpp(1:obj%NxBins) )
  allocate( obj%yLow(1:obj%NyBins) ,obj%yUpp(1:obj%NyBins) )
  allocate( obj%sumw(1:2,1:obj%NxBins,1:obj%NyBins) )
  do ii=1,obj%NxBins
    obj%xLow(ii) = xBins(2*ii-1)
    obj%xUpp(ii) = xBins(2*ii  )
  enddo
  do ii=1,obj%NyBins
    obj%yLow(ii) = yBins(2*ii-1)
    obj%yUpp(ii) = yBins(2*ii  )
  enddo
  obj%Nev = 0
  obj%sumw = 0
  end subroutine


  subroutine collect_2d( obj ,xVal,yVal ,weight )
  class(histo_2d_type) :: obj
  real(kind(1d0)),intent(in) :: xVal,yVal,weight
  integer :: ii,i0,i1,jj,j0,j1
  obj%Nev = obj%Nev+1
  if (weight.eq.0d0) return
  if (xVal.lt.obj%xLow(1)) return
  if (xVal.ge.obj%xUpp(obj%NxBins)) return
  if (yVal.lt.obj%yLow(1)) return
  if (yVal.ge.obj%yUpp(obj%NyBins)) return
  i0 = 1
  i1 = obj%Nxbins
  do while (i0.ne.i1)
    ii = (i0+i1)/2
    if (xVal.lt.obj%xLow(ii)) then
      i1 = ii-1
    elseif (xVal.lt.obj%xUpp(ii)) then
      i1 = ii
      exit
    elseif (xVal.ge.obj%xUpp(ii)) then
      i0 = ii+1
    endif
  enddo
  if (xVal.lt.obj%xLow(i1).or.obj%xUpp(i1).le.xVal) return
  j0 = 1
  j1 = obj%Nybins
  do while (j0.ne.j1)
    jj = (j0+j1)/2
    if (yVal.lt.obj%yLow(jj)) then
      j1 = jj-1
    elseif (yVal.lt.obj%yUpp(jj)) then
      j1 = jj
      exit
    elseif (yVal.ge.obj%yUpp(jj)) then
      j0 = jj+1
    endif
  enddo
  if (yVal.lt.obj%yLow(j1).or.obj%yUpp(j1).le.yVal) return
  obj%sumw(1:2,i1,j1) = obj%sumw(1:2,i1,j1) + [weight,weight*weight]
  end subroutine


  subroutine write_2d( obj ,filename ,frmt ,factor ,writeUnit )
  class(histo_2d_type) :: obj
  character(*),intent(in) :: filename
  character(*),intent(in),optional :: frmt
  real(kind(1d0)),intent(in),optional :: factor
  integer,intent(in),optional :: writeUnit
  integer :: ii,jj,wUnit
  real(kind(2d0)) :: h0,h1,h2,h3,h4
  h0 = obj%Nev
  wUnit=59;if(present(writeUnit))wUnit=writeUnit
  open(wUnit,file=trim(filename),status='new')
  !This is f2008: open(newunit=wUnit,file=trim(filename),status='new')
  select case (frmt)
  case default
    do ii=1,obj%NxBins
    do jj=1,obj%NyBins
      call compute
      write(wUnit,'(6e16.8)') obj%xLow(ii),obj%yLow(jj),obj%xUpp(ii),obj%yUpp(jj),h1,h2
    enddo
    enddo
  case ('gnuplot')
    do ii=1,obj%NxBins
      do jj=1,obj%NyBins
        call compute
        write(wUnit,'(3e16.8)') obj%xLow(ii),obj%yLow(jj),h1
        write(wUnit,'(3e16.8)') obj%xLow(ii),obj%yUpp(jj),h1
      enddo
      write(wUnit,*)
      do jj=1,obj%NyBins
        call compute
        write(wUnit,'(3e16.8)') obj%xUpp(ii),obj%yLow(jj),h1
        write(wUnit,'(3e16.8)') obj%xUpp(ii),obj%yUpp(jj),h1
      enddo
      write(wUnit,*)
    enddo
  end select
  close(wUnit)
  contains
    subroutine compute
    if (h0.le.0) then
      h1 = 0
      h2 = 0
    else
      h1 = obj%sumw(1,ii,jj)/h0
      if (h0.le.1) then
        h2 = abs(h1)
      else
        h2 = (obj%sumw(2,ii,jj)/h0 - h1*h1)/(h0-1)
        if (h2.lt.0) then
          h2 = abs(h1)
        else
          h2 = sqrt(h2)
        endif
      endif
    endif
    h1 = h1/(obj%xUpp(ii)-obj%xLow(ii))/(obj%yUpp(jj)-obj%yLow(jj))
    h2 = h2/(obj%xUpp(ii)-obj%xLow(ii))/(obj%yUpp(jj)-obj%yLow(jj))
    if (present(factor)) then
      h1 = h1*factor
      h2 = h2*factor
    endif
    end subroutine
  end subroutine


  function mass(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt
  rslt = sqrt(abs( (p(0)+p(3))*(p(0)-p(3)) - p(1)*p(1) - p(2)*p(2) ))
  end function

  function pTrans(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt
  rslt = sqrt( p(1)*p(1)+p(2)*p(2) )
  end function

  function ETrans(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt
  rslt = sqrt(abs( (p(0)+p(3))*(p(0)-p(3)) ))
  end function

  function rapidity(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt,aa,bb
  bb = abs(p(0))
  aa = bb+p(3)
  bb = bb-p(3)
  if (aa*bb.gt.0) then
    rslt = log(aa/bb)/2
  else
    rslt = 99
  endif
  if (p(0).lt.0) rslt =-rslt
  end function

  function pseudoRap(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt,pTsq
  pTsq = p(1)*p(1) + p(2)*p(2)
  if (pTsq.gt.0) then
    rslt = sqrt(pTsq+p(3)*p(3)) + abs(p(3))
    rslt = log( rslt*rslt/pTsq )/2
  else
    rslt = 99
  endif
  if (p(3).lt.0) rslt =-rslt
  end function

  function phi(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt,pTsq
  pTsq = p(1)*p(1) + p(2)*p(2)
  if (pTsq.le.0) then
    if (p(3).ge.0) then ;rslt = 0
                   else ;rslt = r1PI
    endif
  else
    rslt = acos(p(1)/sqrt(pTsq))
    if (p(2).lt.0) rslt=r2PI-rslt
  endif
  end function

  function deltaR(p,q) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3),q(0:3)
  real(kind(1d0)) :: rslt
  rslt = abs( phi(p) - phi(q) )
  rslt = min( rslt ,r2PI-rslt )
  rslt = sqrt( rslt*rslt + (rapidity(p)-rapidity(q))**2 )
  end function

  function deltaPhi(p,q) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3),q(0:3)
  real(kind(1d0)) :: rslt
  rslt = abs( phi(p) - phi(q) )
  rslt = min( rslt ,r2PI-rslt )
  end function

  function theta(p) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3)
  real(kind(1d0)) :: rslt
  rslt = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
  if (rslt.le.0) then
    rslt = 0
  else
    rslt = acos( p(3)/sqrt(rslt) )
    rslt = rslt/r1PI*180
  endif
  end function

  function angle(p,q) result(rslt)
  real(kind(1d0)),intent(in) :: p(0:3),q(0:3)
  real(kind(1d0)) :: rslt ,ap,aq
  ap = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
  aq = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
  if (ap.le.0.or.aq.le.0) then
    rslt = 0
  else
    rslt = ( p(1)*q(1) + p(2)*q(2) + p(3)*q(3) )/sqrt(ap*aq)
    rslt = acos(rslt)/r1PI*180
  endif
  end function


  subroutine breit0_init(obj,qq)!,discard)
! Boost momenta to the Breit frame of qq
! qq must have positive energy and negative square
  class(breit0_type) :: obj
  real(kind(1d0)),intent(in) :: qq(0:3)
!  logical,intent(out) :: discard
!  discard = .false.
  obj%sgn = sign(1d0,qq(3))
  obj%rr(0) = qq(1)*qq(1) + qq(2)*qq(2) + qq(3)*qq(3)
  obj%mm = obj%rr(0)-qq(0)*qq(0)
  if (obj%mm.le.0) then
    write(*,*) 'ERROR: Breit frame not defined.'
!    discard = .true.
    obj%sgn = 1
    obj%rr = [1,0,0,0]
    obj%mm = 1
    obj%stsp=0 ; obj%stcp=0 ; obj%ct=1
  else
    obj%mm = sqrt(obj%mm)
    obj%rr(0) = sqrt(obj%rr(0))
    obj%stsp = qq(1)/obj%rr(0) ; obj%stcp = qq(2)/obj%rr(0) ; obj%ct = qq(3)/obj%rr(0)
    obj%rr(1) = obj%stsp*qq(0) ; obj%rr(2) = obj%stcp*qq(0) ; obj%rr(3) = obj%ct*qq(0)
  endif
  obj%st = sqrt(obj%stsp*obj%stsp + obj%stcp*obj%stcp)
  if (obj%st.ne.0) then
    obj%sp = obj%stsp/obj%st
    obj%cp = obj%stcp/obj%st
  else
    obj%sp = 0
    obj%cp = 1
  endif
  obj%ctsp = obj%ct*obj%sp
  obj%ctcp = obj%ct*obj%cp
  end subroutine

  function breit0_act(obj,pp) result(rslt)
  class(breit0_type) :: obj
  real(kind(1d0)),intent(in) :: pp(0:3)
  real(kind(1d0)) :: rslt(0:3),uu(1:3),bb
! boost 
  rslt(0) = ( obj%rr(0)*pp(0) - obj%rr(1)*pp(1) &
            - obj%rr(2)*pp(2) - obj%rr(3)*pp(3) )/obj%mm
  bb = ( pp(0) + rslt(0) )/( obj%rr(0) + obj%mm )
  uu(1:3) = pp(1:3)-bb*obj%rr(1:3)
! rotate
  rslt(1) =   obj%cp*uu(1) -   obj%sp*uu(2)
  rslt(2) = obj%ctsp*uu(1) + obj%ctcp*uu(2) - obj%st*uu(3)
  rslt(3) = obj%stsp*uu(1) + obj%stcp*uu(2) + obj%ct*uu(3)
  rslt(1:3) = obj%sgn*rslt(1:3)
  end function

  function breit0_act_EV(obj,ee,vv) result(pp)
  class(breit0_type) :: obj
  real(kind(1d0)),intent(in) :: ee,vv(1:3)
  real(kind(1d0)) :: pp(0:3)
  pp = obj%act([ee,vv(1),vv(2),vv(3)])
  end function


  subroutine breit_init(obj,qq)!,discard)
! Boost momenta to the Breit frame of qq
! Makes sure that if qq=p1-p2 with p1^2=p2^2=0,
! then the boosted vec{p1},vec{p2} are in the x-z plane.
  class(breit_type) :: obj
  real(kind(1d0)),intent(in) :: qq(0:3)
!  logical,intent(out) :: discard
  real(kind(1d0)) :: aa,bb,cc
!  discard = .false.
  aa = qq(1)*qq(1)+qq(2)*qq(2)
  if (aa.le.0) then
    write(*,*) 'WARNING: Breit frame not defined. a'
    aa = 1
!    discard = .true.
    return
  endif
  aa = sqrt(qq(1)*qq(1)+qq(2)*qq(2))
  obj%sp = qq(1)/aa
  obj%cp = qq(2)/aa 
  aa = obj%sp*qq(1) + obj%cp*qq(2)
  bb = aa*aa + qq(3)*qq(3) - qq(0)*qq(0)
  if (bb.le.0) then
    write(*,*) 'WARNING: Breit frame not defined. b'
    bb = 1
!    discard = .true.
    return
  endif
  cc = qq(0)-qq(3)
  if (cc.eq.0) then
    write(*,*) 'WARNING: Breit frame not defined. c'
    cc = 1
!    discard = .true.
    return
  endif
  bb = sqrt(bb)
  obj%c00 = qq(0)/bb+bb/cc ;obj%c01 =-aa/bb ;obj%c03 =-qq(3)/bb-bb/cc
  obj%c10 =-aa/cc                           ;obj%c13 = aa/cc
  obj%c30 = qq(0)/bb       ;obj%c31 =-aa/bb ;obj%c33 =-qq(3)/bb
  end subroutine

  function breit_act(obj,pp) result(rslt)
  class(breit_type) :: obj
  real(kind(1d0)),intent(in) :: pp(0:3)
  real(kind(1d0)) :: rslt(0:3), p1
  p1 = obj%sp*pp(1) + obj%cp*pp(2)
  rslt(0) = obj%c00*pp(0) + obj%c01*p1 + obj%c03*pp(3)
  rslt(1) = obj%c10*pp(0) +         p1 + obj%c13*pp(3)
  rslt(2) = -obj%cp*pp(1) + obj%sp*pp(2)
  rslt(3) = obj%c30*pp(0) + obj%c31*p1 + obj%c33*pp(3)
  end function

  function breit_act_EV(obj,ee,vv) result(pp)
  class(breit_type) :: obj
  real(kind(1d0)),intent(in) :: ee,vv(1:3)
  real(kind(1d0)) :: pp(0:3)
  pp = obj%act([ee,vv(1),vv(2),vv(3)])
  end function


  subroutine sort_small2big( ll ,aa ,Nsize )
  real(kind(1d0)),intent(inout) :: aa(:)
  integer        ,intent(out  ) :: ll(:)
  integer,intent(in),optional :: Nsize
  integer :: ii,jj,mm,nn
  real(kind(1d0)) :: bb
  if (present(Nsize)) then
    nn = Nsize
  else
    nn = size(aa,1)
  endif
  do ii=1,nn
    ll(ii) = ii
  enddo
  do ii=2,nn
    mm = ll(ii)
    bb = aa(ii)
    jj = ii-1
    do
      if (jj.lt.1) exit
      if (aa(jj).le.bb) exit
      ll(jj+1) = ll(jj)
      aa(jj+1) = aa(jj)
      jj = jj-1
    enddo
    ll(jj+1) = mm
    aa(jj+1) = bb
  enddo
  end subroutine

  subroutine sort_big2small( ll ,aa ,Nsize )
  real(kind(1d0)),intent(inout) :: aa(:)
  integer        ,intent(out  ) :: ll(:)
  integer,intent(in),optional :: Nsize
  integer :: ii,jj,mm,nn
  real(kind(1d0)) :: bb
  if (present(Nsize)) then
    nn = Nsize
  else
    nn = size(aa,1)
  endif
  do ii=1,nn
    ll(ii) = ii
  enddo
  do ii=nn-1,1,-1
    mm = ll(ii)
    bb = aa(ii)
    jj = ii+1
    do
      if (jj.gt.nn) exit
      if (aa(jj).le.bb) exit 
      ll(jj-1) = ll(jj)
      aa(jj-1) = aa(jj)
      jj = jj+1
    enddo
    ll(jj-1) = mm
    aa(jj-1) = bb
  enddo
  end subroutine


end module


