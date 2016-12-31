!  Tri.f90 

!****************************************************************************
!
!  PROGRAM: Tri
!
!  PURPOSE:  Solves tridiagonal linear systems.
!
!****************************************************************************

!
! Numerical Mathematics and Computing, Fifth Edition
! Ward Cheney & David Kincaid
! Brooks/Cole Publ. Co.
! Copyright (c) 2003.  All rights reserved.
! For educational use with the Cheney-Kincaid textbook.
! Absolutely no warranty implied or expressed.
!
! Section 7.3
!
! File: tri.f90
!
! Solves tridiagonal linear systems (tri)
! the solution for each system is x(i)=1

program main
      implicit none
      
interface                                                       
      subroutine tri(n,a,d,c,b,x)    
      implicit none
      integer, intent(in)::n                                    
      real, dimension(n), intent(inout):: a,d,c,b                  
      real, dimension(n), intent(out)::x                        
      end subroutine tri                                        
end interface 

      real, dimension (10):: a,d,c,b,x
      integer ::n = 10   
      integer :: i


      do i=1,n 
      d(i) = 2.0    
      a(i) = 0.5
      c(i) = 0.5      
      b(i) = 3.0 
      end do
      b(1) = 2.5
      b(n) = 2.5     
      call tri(n,a,d,c,b,x) 
      print*, "subroutine tri" 
      print"(f22.14)",(x(i),i=1,n)  
      do  i=1,n 
      d(i) = 2.0    
      c(i) = 0.5   
      b(i) = 3.0 
      end do
      b(1) = 2.5
      b(n) = 2.5     
      call tri(n,c,d,c,b,x) 
      print*, "subroutine tri again" 
      print "(f22.14)",(x(i),i=1,n)  
    end program main

    subroutine tri(n,a,d,c,b,x)    
      implicit none
      integer, intent(in)::n                                       
      real, dimension(n), intent(inout)::a,d,c,b                      
      real, dimension(n), intent(out):: x                          
      integer ::i                                                  
      real :: xmult                                                
      do i = 2,n                                                   
        xmult = a(i-1)/d(i-1)                                      
        d(i) = d(i) - xmult*c(i-1)                                 
        b(i) = b(i) - xmult*b(i-1)                                 
      end do                                                       
      x(n) = b(n)/d(n)                                             
      do i = n-1,1,-1                                              
        x(i) = (b(i) - c(i)*x(i+1))/d(i)                           
      end do                                                       
    end subroutine tri   
  

