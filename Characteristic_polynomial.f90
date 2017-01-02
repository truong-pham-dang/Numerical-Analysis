!  Characteristic_polynomial.f90 
!

!****************************************************************************
!
!  PROGRAM: Characteristic_polynomial
!
!  PURPOSE:  To compute characteristic polynomial of a real 3x3 matrix.
!
!  REFERENCE: charpoly3.pdf, part 1 and section 2.1
!
!****************************************************************************

    program Characteristic_polynomial

    implicit none

    ! Variables
    real :: a(3,3), det
    integer :: i, j 
    real :: c1, c2, c3
    real :: temp(3,3)

    ! Body of Characteristic_polynomial
    
    a(1,:) = (/3.0, 1.0, 1.0/)
    a(2,:) = (/2.0, 4.0, 2.0/)
    a(3,:) = (/1.0, 1.0, 3.0/)
    
    do i=1, 3, 1
     if(i==2) then 
        write(*,*) 'A = [',(a(i,j),j=1,3,1),']'
     else 
        write(*,*) '    [',(a(i,j),j=1,3,1),']'
     end if
    end do
    
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
       + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3))  &
       + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
    print *, '|A| =', det
    
    c1 = 0.0
    do i = 1, 3
        c1 = c1 + a(i,i)
    enddo  
    c1 = -1.0 * c1
    
    temp = 0.0
    temp = matmul(a,a) + c1*a
    
    c2 = 0.0
    do i = 1, 3
        c2 = c2 + temp(i,i)
    enddo
    c2 = (-1.0/2.0)*c2
    
    c3 = (-1)**3 * det
    
    print*,'The characteristic polynomial is:'
    print*,'X^3 + c1*X^2 + c2*X + c3, where'
    print*,'c1 =', c1
    print*,'c2 =', c2
    print*,'c3 =', c3

    end program Characteristic_polynomial

