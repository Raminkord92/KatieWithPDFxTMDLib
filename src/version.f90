module katie_version
  implicit none
  logical,save :: initd=.false.
contains
  subroutine version
  if (initd) return
  initd = .true.
                 !543210987654321098765432109876543210012345678901234567890123456789012345
  write(*,'(A)') '########################################################################'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '#                         You are using KaTie                          #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '#       for event generation with K_T-dependent initial states.        #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
  write(*,'(A)') '#   date: 2021-11-28                                                   #'
  write(*,'(A)') '#    git: http://bitbucket.org/hameren/katie                           #'
  write(*,'(A)') '#                                                                      #'
  write(*,'(A)') '# Please cite                                                          #'
  write(*,'(A)') '#   A. van Hameren, Comput.Phys.Commun. 224 (2018) 371-380             #'
  write(*,'(A)') '# in publications with results obtained with the help of this program. #'
  write(*,'(A)') '########################################################################'
  end subroutine
end module
