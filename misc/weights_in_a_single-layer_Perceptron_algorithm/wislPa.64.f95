! Last updated: 18-Oct-2012
! Udaya Maurya (udaya_cbscients@yahoo.com, telegram: https://t.me/udy11)
! Source: https://github.com/udy11, https://gitlab.com/udy11
! Program (64-bit) to calculate weights in a
! single-layer Perceptron algorithm

! INPUT:
! See the input file provided and make inputs accordingly
! The function f defined in this program is also an input
! Many other important functions are also provided but not used

! Simply keep on pressing enter to see the weights
! in next iterations
! The formula for j^th input is:
! y(j)=f(b+sum_i(w(i)*x(i,j)))
! Then weights are changed as:
! w(i)=w(i)+r*(d(j)-y(j))*x(i,j)

! Sample given is from Wikipedia page of Perceptron
! It uses 4 inputs of NAND as Training Set
! The first entry is held constant at 1

    real*8 x(10,10),d(10),w(10),y(10),b,r
    open(63,file='wislPa_in.txt')
    read(63,*) n,k,m,b,r
    read(63,*) (w(i),i=1,k)
    do j=1,n
       read(63,*) (x(i,j),i=1,k),d(j)
	enddo
	ni=0
47  do l=1,m
    do j=1,n
        y(j)=b+w(1)*x(1,j)
        do i=2,k
            y(j)=y(j)+w(i)*x(i,j)
        enddo
        y(j)=f(y(j))     ! Change the function here
        do i=1,k
           w(i)=w(i)+r*(d(j)-y(j))*x(i,j)
        enddo
	enddo
	write(*,*) (w(i),i=1,k)
	enddo
	ni=ni+m
	write(*,*) "Total iterations made:",ni
    read(*,*)
    goto 47
    end

! Threshold Function
    real*8 function f(z)
    real*8 z,t
    t=0.5d0          ! Change the threshold here
    if(z<t) then
    	f=0.0d0
    else
    	f=1.0d0
	endif
	endfunction
	
! Sigmoid Function
    real*8 function g(z)
    real*8 z,a
    a=1.0d0         ! Slope Parameter
    g=1.0d0/(1.0d0+dexp(-a*z))
	endfunction

! Signum Function
    real*8 function p(z)
    real*8 z
    if(z<0.0d0) then
    	p=-1.0d0
    else
    	p=1.0d0
	endif
	endfunction
	
! Hyperbolic Tangent Function
    real*8 function q(z)
    real*8 z
    q=dtanh(z)
	endfunction
