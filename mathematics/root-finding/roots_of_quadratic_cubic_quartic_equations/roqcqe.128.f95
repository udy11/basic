! Last updated: 08-Mar-2015
! Udaya Maurya (udaya_cbscients@yahoo.com, telegram: https://t.me/udy11)
! Source: https://github.com/udy11, https://gitlab.com/udy11
! Subroutines (128-bit) to solve quadratic, cubic and
! quartic equations with complex coefficients

! ALL YOU NEED TO DO:
! Call qdreqn_128(a) to solve quadratic equation
!   a(1) + a(2) * x + a(3) * x**2 == 0
!   or cbceqn_128(a) to solve cubic equation
!   a(1) + a(2) * x + a(3) * x**2 + a(4) * x**3 == 0
!   or qrteqn_128(a) to solve quartic equation
!   a(1) + a(2) * x + a(3) * x**2 + a(4) * x**3 + a(5) * x**4 == 0
!   where a(i) and x are complex numbers

! NOTE:
! Do not forget to properly declare arrays
! The cbceqn_128 and qrteqn_128 are not working
!   as tested with gfortran in Windows 8.1 x64
!   most probably due to some problem with
!   taking cube roots in gfortran

    complex*32 a,x,b,y,c,z
    dimension a(3),x(2),b(4),y(3),c(5),z(4)

    a(1)=2.q0; a(2)=(3.q0,-7.q0); a(3)=(0.q0,5.q0)
    b(1)=9.q0; b(2)=0.q0; b(3)=(3.q0,7.q0); b(4)=(0.q0,2.q0)
    c(1)=9.q0; c(2)=0.q0; c(3)=(3.q0,7.q0); c(4)=(0.q0,2.q0); c(5)=(-5.q0,-2.q0)

    call qdreqn_128(a,x)
    write(6,*) x
    call cbceqn_128(b,y)
    write(6,*) y
    call qrteqn_128(c,z)
    write(6,*) z
    end

! Subroutine (128-bit) to solve quadratic equation
! a(1) + a(2) * x + a(3) * x**2 == 0
! a and x are complex arrays
    subroutine qdreqn_128(a,x)
    complex*32 a,x,a2,ds
    dimension a(3),x(2)
    if(a(3)==0.q0) then
        if(a(2)/=0.q0) then
            x(1)=-a(1)/a(2); x(2)=x(1)
        endif
    else
        a2=a(2)*0.5q0
        ds=sqrt(a2*a2-a(3)*a(1))
        x(1)=(-a2+ds)/a(3)
        x(2)=(-a2-ds)/a(3)
    endif
    endsubroutine

! Subroutine (128-bit) to solve cubic equation
! a(1) + a(2) * x + a(3) * x**2 + a(4) * x**3 == 0
! a and x are complex arrays
    subroutine cbceqn_128(a,x)
    complex*32 a,x,a2,ds,a24,a33,d0,d1,cd,cc,u1,u2,tt
    dimension a(4),x(3)
    if(a(4)==0.q0) then
        if(a(3)==0.q0) then
            if(a(2)/=0.q0) then
                x(1)=-a(1)/a(2); x(2)=x(1)
            endif
        else
            a2=a(2)*0.5q0
            ds=sqrt(a2*a2-a(3)*a(1))
            x(1)=(-a2+ds)/a(3)
            x(2)=(-a2-ds)/a(3)
        endif
        x(3)=x(2)
    else
        a24=a(2)*a(4)
        a33=a(3)*a(3)
        d0=a33-3.q0*a24
        d1=(2.q0*a33-9.q0*a24)*a(3)+27.q0*a(4)*a(4)*a(1)
        if(d0==(0.q0,0.q0) .and. d1==(0.q0,0.q0)) then
            x=-0.33333333333333333333333333333333333q0*a(3)/a(4)
            return
        endif
        cd=sqrt(d1*d1-4.q0*d0*d0*d0)
        if(abs(d1+cd)>abs(d1-cd)) then
            cc=((d1+cd)*0.5q0)**0.33333333333333333333333333333333333q0
        else
            cc=((d1-cd)*0.5q0)**0.33333333333333333333333333333333333q0
        endif
        u1=(-0.5q0,0.86602540378443864676372317075293618q0)*cc
        u2=(-0.5q0,-0.86602540378443864676372317075293618q0)*cc
        tt=-0.33333333333333333333333333333333333q0/a(4)
        x(1)=tt*(a(3)+cc+d0/cc)
        x(2)=tt*(a(3)+u1+d0/u1)
        x(3)=tt*(a(3)+u2+d0/u2)
    endif
    endsubroutine

! Subroutine (128-bit) to solve quartic equation
! a(1) + a(2) * x + a(3) * x**2 + a(4) * x**3 + a(5) * x**4 == 0
! a and x are complex arrays
    subroutine qrteqn_128(a,x)
    complex*32 a,x,a2,ds,a24,a33,d0,d1,cd,cc,u1,u2,tt
    complex*32 a44,a35,a55,p,q,qds,q1,s1,qq,ss,b45,sp,qs
    dimension a(5),x(4)
    if(a(5)==0.q0) then
        if(a(4)==0.q0) then
            if(a(3)==0.q0) then
                if(a(2)/=0.q0) then
                    x(1)=-a(1)/a(2); x(2)=x(1)
                endif
            else
                a2=a(2)*0.5q0
                ds=sqrt(a2*a2-a(3)*a(1))
                x(1)=(-a2+ds)/a(3)
                x(2)=(-a2-ds)/a(3)
            endif
            x(3)=x(2)
        else
            a24=a(2)*a(4)
            a33=a(3)*a(3)
            d0=a33-3.q0*a24
            d1=(2.q0*a33-9.q0*a24)*a(3)+27.q0*a(4)*a(4)*a(1)
            if(d0==(0.q0,0.q0) .and. d1==(0.q0,0.q0)) then
                x=-0.33333333333333333333333333333333333q0*a(3)/a(4)
                return
            endif
            cd=sqrt(d1*d1-4.q0*d0*d0*d0)
            if(abs(d1+cd)>abs(d1-cd)) then
                cc=((d1+cd)*0.5q0)**0.33333333333333333333333333333333333q0
            else
                cc=((d1-cd)*0.5q0)**0.33333333333333333333333333333333333q0
            endif
            u1=(-0.5q0,0.86602540378443864676372317075293618q0)*cc
            u2=(-0.5q0,-0.86602540378443864676372317075293618q0)*cc
            tt=-0.33333333333333333333333333333333333q0/a(4)
            x(1)=tt*(a(3)+cc+d0/cc)
            x(2)=tt*(a(3)+u1+d0/u1)
            x(3)=tt*(a(3)+u2+d0/u2)
        endif
    else
        a44=a(4)*a(4)
        a35=a(3)*a(5)
        a55=a(5)*a(5)
        p=(a35-0.375q0*a44)/a55
        q=((0.125q0*a44-0.5q0*a35)*a(4)+a55*a(2))/(a(5)*a55)
        d0=a(3)*a(3)-3.q0*(a(2)*a(4)-4.q0*a(1)*a(5))
        d1=(a(3)*a(3)-4.5q0*(a(4)*a(2)+8.q0*a(1)*a(5)))*a(3)+13.5q0*(a(2)*a(2)*a(5)+a(1)*a44)
        if(d1==(0.q0,0.q0) .and. d0==(0.q0,0.q0)) then
            if(p==(0.q0,0.q0)) then
                x=-0.25q0*a(4)/a(5)
                return
            endif
            ss=0.28867513459481288225457439025097873q0*sqrt(-2.q0*p)
        else
            qds=sqrt(d1*d1-d0*d0*d0)
            if(abs(d1+qds)>abs(d1-qds)) then
                qq=(d1+qds)**0.33333333333333333333333333333333333q0
            else
                qq=(d1-qds)**0.33333333333333333333333333333333333q0
            endif
            ss=0.28867513459481288225457439025097873q0*sqrt(-2.q0*p+(qq+d0/qq)/a(5))
            q1=qq*(-0.5q0,0.86602540378443864676372317075293618q0)
            s1=0.28867513459481288225457439025097873q0*sqrt(-2.q0*p+(q1+d0/q1)/a(5))
            if(abs(s1)>abs(ss)) ss=s1
            q1=qq*(-0.5q0,-0.86602540378443864676372317075293618q0)
            s1=0.28867513459481288225457439025097873q0*sqrt(-2.q0*p+(q1+d0/q1)/a(5))
            if(abs(s1)>abs(ss)) ss=s1
        endif
        b45=-0.25q0*a(4)/a(5)
        sp=-ss*ss-0.5q0*p
        qs=0.25q0*q/ss
        x(1)=b45-ss-sqrt(sp+qs)
        x(2)=b45-ss+sqrt(sp+qs)
        x(3)=b45+ss-sqrt(sp-qs)
        x(4)=b45+ss+sqrt(sp-qs)
    endif
    endsubroutine
