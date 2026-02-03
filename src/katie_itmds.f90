module katie_itmds
  use avh_trivinp
  use avh_lagrangian
  implicit none
  private
  public :: gg1,gg2,gg3,gg4,gg5,gg6,gg7,qg1,qg2,qg3
  public :: pdfs_init,set_tmdpath,add_tmdpdf,alphasFunc,tmdfunc,pdfunc_lhapdf

  integer,parameter :: qg1=1,qg2=2,qg3=3
  integer,parameter :: gg1=4,gg2=5,gg3=6,gg4=7,gg5=8,gg6=9,gg7=10

  type(trivinp_type),save :: tmdMesh(10)
  logical :: loaded(10)=.false.

  integer,save :: lgn2pdf(nulParticle:maxParticle)
  integer,parameter :: gl=0,uQ=1,dQ=2,sQ=3,cQ=4,bQ=5,tQ=6
  integer,parameter :: pdf2lha(-6:6)=[-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6]

  character(256),save :: path

  interface tmdfunc
    module procedure tmdfunc_1,tmdfunc_2
  end interface

contains

  subroutine pdfs_init( Nflavors ,lhaSet  )
  integer,intent(in) :: Nflavors !dummy for standard function format
  character(*),intent(in) :: lhaSet
  lgn2pdf = 1001
  lgn2pdf(gluon )=gl
  lgn2pdf(uQuark)=uQ; lgn2pdf(-uQuark)=-lgn2pdf(uQuark);
  lgn2pdf(dQuark)=dQ; lgn2pdf(-dQuark)=-lgn2pdf(dQuark);
  lgn2pdf(cQuark)=cQ; lgn2pdf(-cQuark)=-lgn2pdf(cQuark);
  lgn2pdf(sQuark)=sQ; lgn2pdf(-sQuark)=-lgn2pdf(sQuark);
  lgn2pdf(tQuark)=tQ; lgn2pdf(-tQuark)=-lgn2pdf(tQuark);
  lgn2pdf(bQuark)=bQ; lgn2pdf(-bQuark)=-lgn2pdf(bQuark);
  call InitPDFsetByName(trim(lhaSet)//'.LHgrid')
  end subroutine

  subroutine set_tmdpath( pathVal )
  character(*),intent(in) :: pathVal
  path = pathVal
  end subroutine

  subroutine add_tmdpdf( label ,fileName )
  character(*),intent(in) :: label,fileName
  select case (label)
  case ('gg1') ;call load(gg1)
  case ('gg2') ;call load(gg2)
  case ('gg3') ;call load(gg3)
  case ('gg4') ;call load(gg4)
  case ('gg5') ;call load(gg5)
  case ('gg6') ;call load(gg6)
  case ('gg7') ;call load(gg7)
  case ('qg1') ;call load(qg1)
  case ('qg2') ;call load(qg2)
  case ('qg3') ;call load(qg3)
  case default
    write(*,*) 'ERROR in itmdf_mod add_tmdpdf: label beyond bounds tmdMesh'
    stop
  end select
  contains
    subroutine load(ii)
    integer,intent(in) :: ii
    call tmdMesh(ii)%read_file(trim(path)//trim(fileName))
    loaded(ii) = .true.
    end subroutine
  end subroutine

  function alphasFunc(xx) result(rslt)
  !(realknd2!),intent(in) :: xx
  !(realknd2!) :: rslt,alphasPDF
  rslt = alphasPDF(xx)
  end function

  function tmdfunc_1( xx ,mu ,kTsq ) result(rslt)
  !(realknd2!),intent(in) :: xx,mu,kTsq
  !(realknd2!) :: rslt(10),logx,logkTsq,logmuSq
  integer :: ii
  rslt = 0
  logx = log(xx)
  logkTsq = log(kTsq)
  logmuSq = log(mu*mu)
  do ii=1,10
    if (loaded(ii)) then
      select case (tmdMesh(ii)%Ndim)
      case default
        rslt(ii) = tmdMesh(ii)%evaluate( logx ,logkTsq ,logmuSq )
        rslt(ii) = rslt(ii)/xx
      case (2)
        rslt(ii) = tmdMesh(ii)%evaluate( logx ,logkTsq )
        rslt(ii) = rslt(ii)/xx
      end select
    endif
  enddo
  end function

  function tmdfunc_2( list ,xx ,mu ,kTsq ) result(rslt)
  integer,intent(in) :: list(:)
  !(realknd2!),intent(in) :: xx,mu,kTsq
  !(realknd2!) :: rslt(10),logx,logkTsq,logmuSq
  integer :: ii
  rslt = 0
  logx = log(xx)
  logkTsq = log(kTsq)
  logmuSq = log(mu*mu)
  do ii=1,size(list)
    select case (tmdMesh(list(ii))%Ndim)
    case default
      rslt(list(ii)) = tmdMesh(list(ii))%evaluate( logx ,logkTsq ,logmuSq )
      rslt(list(ii)) = rslt(list(ii))/xx
    case (2)
      rslt(list(ii)) = tmdMesh(list(ii))%evaluate( logx ,logkTsq )
      rslt(list(ii)) = rslt(list(ii))/xx
    end select
  enddo
  end function

  function pdfunc_lhapdf( iParton ,xx ,mu ) result(rslt)
  integer,intent(in) :: iParton
  !(realknd2!),intent(in) :: xx,mu
  !(realknd2!) :: rslt,list(-6:6)
  rslt = 1
  if (lgn2pdf(iParton).eq.1001) return
  rslt = 0
  call evolvePDF(xx,mu,list)
  rslt = list(pdf2lha(lgn2pdf(iParton)))/xx
  end function

end module


