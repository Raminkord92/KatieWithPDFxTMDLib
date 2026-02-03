module katie_events
  use avh_iounits
  use avh_mathcnst
  implicit none
  private
  public :: init_events ,update_events ,write_accumulated ,write_result

  integer,parameter :: NstpMessage=20

  integer,save :: wUnit=21,mUnit=6,Nmessage=NstpMessage**2,Nevents,Nbuffer
  logical,save :: initPhase=.true.
  !(realknd2!),save :: sumW0tot,sumWevent
  !(realknd2!),save :: sumW0pos,sumW1pos,sumW2pos
  !(realknd2!),save :: sumW0neg,sumW1neg,sumW2neg
  !(realknd2!),save :: targetPrecision,NevParameter,Wthreshold
  !(realknd2!),allocatable,save :: wBuffer(:)
  interface
    function rng_interface() result(rslt)
    !(realknd2!) :: rslt
    end function
  end interface

  procedure(rng_interface),pointer :: rng => NULL()

contains

  
  subroutine init_events( rngenerator ,writeUnit ,messageUnit &
                         ,desiredPrecision ,desiredNevents ,bufferSize )
  procedure(rng_interface) :: rngenerator
  integer,intent(in),optional :: writeUnit,messageUnit,bufferSize
  !(realknd2!),intent(in),optional :: desiredPrecision,desiredNevents
  if (present(writeUnit  )) wUnit = writeUnit
  if (present(messageUnit)) mUnit = messageUnit
  rng => rngenerator
  Nevents = 0
  sumW0tot = 0
  sumWevent = 0
  sumW0pos=0;sumW1pos=0;sumW2pos=0
  sumW0neg=0;sumW1neg=0;sumW2neg=0
  initPhase = .true.
  targetPrecision = 1d-3
    if (present(desiredPrecision)) targetPrecision = desiredPrecision
  NevParameter = 1d5 * targetPrecision**2
    if (present(desiredNevents)) NevParameter = desiredNevents*targetPrecision**2
  if (allocated(wBuffer)) deallocate(wBuffer)
  Nbuffer = 10*1000
    if (present(bufferSize)) Nbuffer = bufferSize
  allocate(wBuffer(0:Nbuffer))
  wBuffer = 0
  end subroutine


  subroutine update_events( discard ,accept ,stopMC ,weight )
  logical,intent(out) :: discard,accept,stopMC
  !(realknd2!),intent(in) :: weight
  !(realknd2!) :: Wevent,absWeight,currentPrecision
  !(realknd2!) :: avePos,varPos,aveNeg,varNeg
  integer :: nn,ii
!
  discard = .false.
  accept = .false.
  stopMC = .false.
!
  if (weight.gt.0) then
    sumW0tot = sumW0tot + 1
    absWeight = weight
    sumW0pos = sumW0pos + 1
    sumW1pos = sumW1pos + absWeight
    sumW2pos = sumW2pos + absWeight*absWeight
  elseif (weight.lt.0) then
    sumW0tot = sumW0tot + 1
    absWeight =-weight
    sumW0neg = sumW0neg + 1
    sumW1neg = sumW1neg + absWeight
    sumW2neg = sumW2neg + absWeight*absWeight
  elseif (weight.eq.0) then
    sumW0tot = sumW0tot + 1
    return
  else
    if (mUnit.ge.0) write(mUnit,*) 'ERROR in events: ' &
      ,'weight=',weight,', discarding event'
    discard = .true.
    return
  endif
!
  call calc_result( avePos ,varPos ,sumW0tot ,sumW1pos ,sumW2pos )
  call calc_result( aveNeg ,varNeg ,sumW0tot ,sumW1neg ,sumW2neg )
  currentPrecision = abs(avePos-aveNeg)
  if (currentPrecision.le.rZRO) then
    currentPrecision = rONE
  else
    currentPrecision = sqrt(varPos+varNeg)/currentPrecision
  endif
!
  if (sumW0pos+sumW0neg.eq.Nmessage) then
    Nmessage = (sqrt(rONE*Nmessage)+NstpMessage)**2
    call write_result
  endif
!
  if (initPhase) then
    nn = sumW0pos+sumW0neg
    wBuffer(nn) = absWeight
    if (nn.eq.Nbuffer) then
      initPhase = .false.
      call sort_small( wBuffer )
      Wthreshold = find_threshold( wBuffer ,NevParameter/(nn*currentPrecision**2) )
      if (mUnit.gt.0) then
        write(mUnit,'(a50,e15.8)') &
          ' MESSAGE: partially un-weighting with Wthreshold =',Wthreshold
        write(mUnit,*)
      endif
      !do ii=1,Nbuffer                                                       !DEBUG
      !  write(88,*) wBuffer(ii-1),wBuffer(ii),1/(wBuffer(ii)-wBuffer(ii-1)) !DEBUG
      !enddo                                                                 !DEBUG
      deallocate(wBuffer)
    endif
  else
    accept = (Wthreshold.le.absWeight)
    if (accept) then
      Wevent = weight
    else
      Wevent = sign(Wthreshold,weight)
      accept = (rng()*Wthreshold.le.absWeight)
    endif
    if (accept) then
      Nevents = Nevents+1
      sumWevent = sumWevent+Wevent
      write(wUnit,'(a13,e24.16)') 'EVENT WEIGHT:',Wevent
    endif
  endif
!
  stopMC = (currentPrecision.le.targetPrecision)
!
  end subroutine


  subroutine write_accumulated
  write(wUnit,'(a17,9e16.9)') 'ACCUMULATED SUMS:' &
       ,rONE*Nevents,sumWevent,sumW0tot &
       ,sumW0pos,sumW1pos,sumW2pos,sumW0neg,sumW1neg,sumW2neg
  end subroutine


  subroutine write_result
  if (mUnit.lt.0) return
  write(mUnit,'(A)') print_whole(sumW0tot) &
                   //print_result(8,sumW0tot,sumW0pos,sumW1pos,sumW2pos)
  if (sumW0neg.ne.0) &
  write(mUnit,'(A)') '   negative:' &
                   //print_result(8,sumW0tot,sumW0neg,sumW1neg,sumW2neg)
  write(mUnit,'(A)')
  end subroutine


  function print_result( nDec ,sum0tot ,sum0 ,sum1 ,sum2 ) result(line)
  !(realknd2!),intent(in) :: sum0tot,sum0,sum1,sum2
  integer,intent(in) :: nDec
  character(34+2*nDec) :: line
  integer :: i1,i2
  !(realknd2!) :: ave,sig,rr,tt,w0
  character(33+2*nDec) :: aa,bb,cc,ee
  call calc_result( ave ,sig ,sum0tot ,sum1 ,sum2 )
  sig = sqrt(sig)
  rr=0 ;if (ave.ne.rZRO) rr=sig/ave
  write(aa,'(e23.16)') ave
  ee(1:2) = '1.'
  ee(3:6) = aa(20:23)
  read(ee(1:6),*) tt
  write(bb,'(f18.16)') abs(sig)/tt
  write(cc,'(f16.12)') abs(rr*100)
  line = print_whole(sum0)//'  ('//aa(3:3+nDec)//'+/-'//bb(2:2+nDec)//')'//ee(3:6) &
       //' '//cc(1:7)//'%'
  end function

  function print_whole( xx ) result(rslt)
  !(realknd2!),intent(in) :: xx
  character(24) :: aa
  character(12) :: rslt
  write(aa,'(f13.0)') xx
  if (aa(9:9).ne.' ') aa(1:9) = aa(2:9)//','
  if (aa(5:5).ne.' ') aa(1:5) = aa(2:5)//','
  rslt = aa(1:12)
  end function


  subroutine calc_result( mean ,var ,sum0 ,sum1 ,sum2 )
  !(realknd2!),intent(out) :: mean ,var
  !(realknd2!),intent(in) :: sum0,sum1,sum2
  !(realknd2!) :: hh
  if (sum0.le.rZRO.or.sum1.le.rZRO) then
    mean = 0
    var = 0
  elseif (sum0.le.rONE) then
    mean = sum1
    var = mean*mean
  else
    mean = sum1/sum0
    hh  = mean*mean
    var = sum2/sum0 - hh
    if (var.le.rZRO) then
      var = hh
    else
      var = var/(sum0-1)
    endif
  endif
  end subroutine


  function find_threshold(ww,rho) result(rslt)
! Solves for  w  from the equation
! /w                  /wMax             /wMax
! | (x/w)*h(x) dx  +  | h(x) dx  =  rho | h(x) dx
! /0                  /w                /0
! where h(x) = sum( theta(w_{i-1}<x<w_i)/(w_i-w_{i-1}) )
!
  !(realknd2!),intent(in) :: ww(0:),rho
  !(realknd2!) :: rslt ,aa,bb,cc,xx
  integer :: i0,i1,ii,nn
  nn = ubound(ww,1)
  if (rho.ge.rONE) then
    rslt = rZRO
  elseif (rho.le.rZRO) then
    rslt = ww(nn)
  else
    i0 = 0
    i1 = nn
    do while (i1-i0.gt.1)
      ii = (i0+i1)/2
      xx = ( sum(ww(1:ii-1))/ww(ii)+rHLF + (nn-ii) )/nn
      if (rho.ge.xx) then
        i1 = ii
      else
        i0 = ii
      endif
    enddo
    aa = 1/(ww(i1)-ww(i0))
    bb = rHLF + ww(i1)*aa + ((1-rho)*nn-i1)
    cc = sum(ww(1:i0))
    rslt = (bb+sqrt(bb*bb+4*aa*cc))/(2*aa)
  endif
  end function


  subroutine sort_small( xx )
  !(realknd2!),intent(inout) :: xx(:)
  integer :: ll,ii,jj,ir,nn
  !(realknd2!) :: xxb
  nn = size(xx)
  if (nn.le.1) return
  ll = nn/2 + 1
  ir = nn
  do
    if (ll.gt.1) then 
      ll = ll - 1
      xxb = xx(ll)
    else
      xxb = xx(ir)
      xx(ir) = xx(1)
      ir = ir - 1
      if (ir.eq.1) then
        xx(1) = xxb
        exit
      endif
    endif
    ii = ll
    jj = ll + ll
    do while (jj.le.ir)
      if (jj.lt.ir) then
        if (xx(jj).lt.xx(jj+1)) jj = jj + 1
!        if (xx(jj).gt.xx(jj+1)) jj = jj + 1
      endif
      if (xxb.lt.xx(jj)) then
!      if (xxb.gt.xx(jj)) then
        xx(ii) = xx(jj)
        ii = jj
        jj = jj + jj
      else
        jj = ir + 1
      endif
    enddo
    xx(ii) = xxb
  enddo
  end subroutine


end module


