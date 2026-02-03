module katie_partlumi
  implicit none
  private

! Before using the functions, the number of flavors needs to be set
!   call set_Nflavors( Nflavors )
  public :: set_Nflavors

! All functions combine evaluated values of pdfs, and have format
!
!   function pdf_A_B( pdfA,pdfB ) result(rslt)
!   real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
!   real(realknd2) :: rslt
!
! quarks are assumed to be enumerated in the mass-ordered convention
! t~=-6 ,b~=-5 ,c~=-4 ,s~=-3 ,d~=-2 ,u~=-1 ,g=0 ,u=1 ,d=2 ,s=3 ,c=4 ,b=5 ,t=6

! Combining pdfs assuming (anti) quarks are equivalent (no EW interaction)
  public :: pdf_g_g ,pdf_g_q ,pdf_q_g ,pdf_q_q ,pdf_q_qb ,pdf_q_r ,pdf_q_rb

! For example, pdf_q_r returns for Nflavors=4
!    u(A)*(        d(B) + s(B) + c(B) ) + u~(A)*(         d~(B) + s~(B) + c~(B) )
!  + d(A)*( u(B)        + s(B) + c(B) ) + d~(A)*( u~(B)         + s~(B) + c~(B) )
!  + s(A)*( u(B) + d(B)        + c(B) ) + d~(A)*( u~(B) + d~(B)         + c~(B) )
!  + c(A)*( u(B) + d(B) + s(B)        ) + d~(A)*( u~(B) + d~(B) + s~(B)         )
!
! For example, pdf_q_rb returns for Nflavors=4
!    u(A)*(         d~(B) + s~(B) + c~(B) ) + u~(A)*(        d(B) + s(B) + c(B) )
!  + d(A)*( u~(B)         + s~(B) + c~(B) ) + d~(A)*( u(B)        + s(B) + c(B) )
!  + s(A)*( u~(B) + d~(B)         + c~(B) ) + d~(A)*( u(B) + d(B)        + c(B) )
!  + c(A)*( u~(B) + d~(B) + s~(B)         ) + d~(A)*( u(B) + d(B) + s(B)        )
!
! For example, pdf_q_qb returns for Nflavors=3
!   u(A)*u~(B) + u~(A)*u(B)  +  d(A)*d~(B) + d~(A)*d(B)  +  s(A)*s~(B) + s~(A)*s(B)

! Combining pdfs taking into account inequivalence of (anti) up/down quarks
  public :: pdf_g_u  ,pdf_g_ub  ,pdf_g_d   ,pdf_g_db  &
           ,pdf_u_g  ,pdf_ub_g  ,pdf_d_g   ,pdf_db_g  &
           ,pdf_u_u  ,pdf_ub_ub ,pdf_d_d   ,pdf_db_db &
           ,pdf_u_ub ,pdf_ub_u  ,pdf_d_db  ,pdf_db_d  &
           ,pdf_u_d  ,pdf_d_u   ,pdf_u_db  ,pdf_db_u  &
           ,pdf_ub_d ,pdf_d_ub  ,pdf_ub_db ,pdf_db_ub

! For example, pdf_u_db returns for Nflavors=4
!   ( u(A)+c(A) )*( d~(B) + s~(B) )
! For example, pdf_u_db returns for Nflavors=3
!   ( u(A) )*( d~(B) + s~(B) )

! up-type and down-type quarks in the mass-ordered (LHAPDF) convention
  integer,parameter :: uTyp(-3:3)=[-6,-4,-1,0,1,4,6]
  integer,parameter :: dTyp(-3:3)=[-5,-3,-2,0,2,3,5]

! number of up/down-type quarks as function of Nflavors
  integer,parameter :: NuTyp(0:6)=[0,1,1,1,2,2,3]
  integer,parameter :: NdTyp(0:6)=[0,0,1,2,2,3,3]

  integer,parameter :: realknd2=kind(1d0)

  integer,save :: Nflavors

contains


  subroutine set_Nflavors( Nflavors_in )
  integer,intent(in) :: Nflavors_in
  Nflavors = Nflavors_in
  end subroutine
  

  function pdf_g_g( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  rslt = pdfA(0)*pdfB(0)
  end function

  function pdf_g_q( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfB(uTyp(ii)) + pdfB(uTyp(-ii))
  enddo
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfB(dTyp(ii)) + pdfB(dTyp(-ii))
  enddo
  rslt = pdfA(0)*rslt
  end function

  function pdf_g_u( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfB(uTyp(ii))
  enddo
  rslt = pdfA(0)*rslt
  end function

  function pdf_g_ub( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfB(uTyp(-ii))
  enddo
  rslt = pdfA(0)*rslt
  end function

  function pdf_g_d( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfB(dTyp(ii))
  enddo
  rslt = pdfA(0)*rslt
  end function

  function pdf_g_db( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfB(dTyp(-ii))
  enddo
  rslt = pdfA(0)*rslt
  end function

  function pdf_q_g( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(ii)) + pdfA(uTyp(-ii))
  enddo
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(ii)) + pdfA(dTyp(-ii))
  enddo
  rslt = rslt*pdfB(0)
  end function

  function pdf_u_g( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(ii))
  enddo
  rslt = rslt*pdfB(0)
  end function

  function pdf_ub_g( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(-ii))
  enddo
  rslt = rslt*pdfB(0)
  end function

  function pdf_d_g( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(ii))
  enddo
  rslt = rslt*pdfB(0)
  end function

  function pdf_db_g( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(-ii))
  enddo
  rslt = rslt*pdfB(0)
  end function

  function pdf_q_q( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(ii))*pdfB(uTyp(ii)) &
                + pdfA(uTyp(-ii))*pdfB(uTyp(-ii))
  enddo
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(ii))*pdfB(dTyp(ii)) &
                + pdfA(dTyp(-ii))*pdfB(dTyp(-ii))
  enddo
  end function

  function pdf_u_u( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(ii))*pdfB(uTyp(ii))
  enddo
  end function

  function pdf_ub_ub( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(-ii))*pdfB(uTyp(-ii))
  enddo
  end function

  function pdf_d_d( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(ii))*pdfB(dTyp(ii))
  enddo
  end function

  function pdf_db_db( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(-ii))*pdfB(dTyp(-ii))
  enddo
  end function

  function pdf_q_qb( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(ii))*pdfB(uTyp(-ii)) &
                + pdfA(uTyp(-ii))*pdfB(uTyp(ii))
  enddo
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(ii))*pdfB(dTyp(-ii)) &
                + pdfA(dTyp(-ii))*pdfB(dTyp(ii))
  enddo
  end function

  function pdf_u_ub( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(ii))*pdfB(uTyp(-ii))
  enddo
  end function

  function pdf_ub_u( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    rslt = rslt + pdfA(uTyp(-ii))*pdfB(uTyp(ii))
  enddo
  end function

  function pdf_d_db( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(ii))*pdfB(dTyp(-ii))
  enddo
  end function

  function pdf_db_d( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt
  integer :: ii
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    rslt = rslt + pdfA(dTyp(-ii))*pdfB(dTyp(ii))
  enddo
  end function

  function pdf_q_r( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh,gg
  integer :: ii,jj
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    hh = 0
    gg = 0
    do jj=1,NuTyp(Nflavors); if(jj.eq.ii)cycle
      hh = hh + pdfB(uTyp( jj))
      gg = gg + pdfB(uTyp(-jj))
    enddo
    do jj=1,NdTyp(Nflavors)
      hh = hh + pdfB(dTyp( jj))
      gg = gg + pdfB(dTyp(-jj))
    enddo
    rslt = rslt +  pdfA(uTyp(ii))*hh + pdfA(uTyp(-ii))*gg
  enddo
  do ii=1,NdTyp(Nflavors)
    hh = 0
    gg = 0
    do jj=1,NuTyp(Nflavors)
      hh = hh + pdfB(uTyp( jj))
      gg = gg + pdfB(uTyp(-jj))
    enddo
    do jj=1,NdTyp(Nflavors); if(jj.eq.ii)cycle
      hh = hh + pdfB(dTyp( jj))
      gg = gg + pdfB(dTyp(-jj))
    enddo
    rslt = rslt + pdfA(dTyp(ii))*hh + pdfA(dTyp(-ii))*gg
  enddo
  end function

  function pdf_q_rb( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh,gg
  integer :: ii,jj
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    hh = 0
    gg = 0
    do jj=1,NuTyp(Nflavors); if(jj.eq.ii)cycle
      hh = hh + pdfB(uTyp( jj))
      gg = gg + pdfB(uTyp(-jj))
    enddo
    do jj=1,NdTyp(Nflavors)
      hh = hh + pdfB(dTyp( jj))
      gg = gg + pdfB(dTyp(-jj))
    enddo
    rslt = rslt +  pdfA(uTyp(ii))*gg + pdfA(uTyp(-ii))*hh
  enddo
  do ii=1,NdTyp(Nflavors)
    hh = 0
    gg = 0
    do jj=1,NuTyp(Nflavors)
      hh = hh + pdfB(uTyp( jj))
      gg = gg + pdfB(uTyp(-jj))
    enddo
    do jj=1,NdTyp(Nflavors); if(jj.eq.ii)cycle
      hh = hh + pdfB(dTyp( jj))
      gg = gg + pdfB(dTyp(-jj))
    enddo
    rslt = rslt + pdfA(dTyp(ii))*gg + pdfA(dTyp(-ii))*hh
  enddo
  end function

  function pdf_u_d( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    hh = 0
    do jj=1,NdTyp(Nflavors)
      hh = hh + pdfB(dTyp(jj))
    enddo
    rslt = rslt + pdfA(uTyp(ii))*hh
  enddo
  end function

  function pdf_d_u( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    hh = 0
    do jj=1,NuTyp(Nflavors)
      hh = hh + pdfB(uTyp(jj))
    enddo
    rslt = rslt + pdfA(dTyp(ii))*hh
  enddo
  end function

  function pdf_u_db( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    hh = 0
    do jj=1,NdTyp(Nflavors)
      hh = hh + pdfB(dTyp(-jj))
    enddo
    rslt = rslt + pdfA(uTyp(ii))*hh
  enddo
  end function

  function pdf_db_u( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    hh = 0
    do jj=1,NuTyp(Nflavors)
      hh = hh + pdfB(uTyp(jj))
    enddo
    rslt = rslt + pdfA(dTyp(-ii))*hh
  enddo
  end function

  function pdf_ub_d( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    hh = 0
    do jj=1,NdTyp(Nflavors)
      hh = hh + pdfB(dTyp(jj))
    enddo
    rslt = rslt + pdfA(uTyp(-ii))*hh
  enddo
  end function

  function pdf_d_ub( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    hh = 0
    do jj=1,NuTyp(Nflavors)
      hh = hh + pdfB(uTyp(-jj))
    enddo
    rslt = rslt + pdfA(dTyp(ii))*hh
  enddo
  end function

  function pdf_ub_db( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NuTyp(Nflavors)
    hh = 0
    do jj=1,NdTyp(Nflavors)
      hh = hh + pdfB(dTyp(-jj))
    enddo
    rslt = rslt + pdfA(uTyp(-ii))*hh
  enddo
  end function

  function pdf_db_ub( pdfA,pdfB ) result(rslt)
  real(realknd2),intent(in) :: pdfA(-6:6),pdfB(-6:6)
  real(realknd2) :: rslt ,hh
  integer :: ii,jj
  rslt = 0
  do ii=1,NdTyp(Nflavors)
    hh = 0
    do jj=1,NuTyp(Nflavors)
      hh = hh + pdfB(uTyp(-jj))
    enddo
    rslt = rslt + pdfA(dTyp(-ii))*hh
  enddo
  end function


end module


