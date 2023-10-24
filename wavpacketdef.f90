!Function to generate the wavefunction values
function psi(x) result(wave)
    real, intent(in)::x
    real::mu,a0,a,V0,x0,sigma0,k0
    complex::wave
    
    mu=1
    a0=0.5
    a=1
    V0=50
    x0=5
    sigma0=1
    k0=10
    n=1

    wave=n*EXP(-(x-x0)**2/(2*sigma0**2))*EXP(-COMPLEX(0,1)*k0*x)
end function

!Function to define the potential at different points in space
function potential(x) result(pot)
    real, intent(in)::x
    real::mu,a0,a,V0,x0,sigma0,k0,pot
    
    mu=1
    a0=0.5
    a=1
    V0=50
    x0=5
    sigma0=1
    k0=10
    n=1

    pot=-V0/(1+EXP((x-a)/a0))+(hbar**2*l*(l+1))/(2*mu*x**2)
end function
program main
    implicit none

    !This part calculates the constant position space
    integer, parameter :: dp = kind(0.0d0)
    real(dp)::nump
    real::hbar,potential,range
    real::delx,dt,delta,l,Hmax,tau
    integer::j,bigj,ii,npoly,num,n,jj
    real, allocatable :: space(:),pot(:,:),kinetic(:,:),H(:,:)
    complex, allocatable :: wave(:),jntau(:),psit(:,:),psishebyexp(:,:),interm(:)
    complex::psi
    
    
    num=200
    nump=200.0
    hbar=1
    range=20.0
    delx=range/num
    dt=0.01
    delta=.1
    l=2.0
    bigj=num
    Hmax = 0+1/(delx**2)
    tau=Hmax*dt

    allocate(space(num))
    do j=0,num-1
        space(j)=(j-num/2)*delx
    end do

    !Now we move on to look at how to construct out wavefunction
    allocate(wave(num))
    do j=0,num-1
        wave(j)=psi(space(j))
    end do

    !Now with this we can begin to generate the chebyshev polynomials with the potential and kinetic energy
    allocate(jntau(14))
    do j=0, 13
        jntau(j)=2*COMPLEX(0,-1)**j*BESSEL_JN(j,tau)
    end do
    npoly=SIZE(jntau)

    !Now we look to generate both our potential and Kinetic energy matrices
    !We use here the 3 point method since it got the best estimation of the conservation of mean energy
    allocate(pot(num,num))
    do j=0,num-1
        if(0<space(j))then
            pot(j,j)=potential(space(j))
        end if
    end do

    !Now we look to define and creat the kinetic energy, here we don't need a function to create the matrix
    allocate(kinetic(num,num))
    do j=0,num-1
        do ii=0,num-1
            if(ii==j) then
                kinetic(ii,j)=1/(delx**2)
            else if(ii-j==ABS(1)) then
                kinetic(ii,j)=-1/(2*delx**2)
            end if
        end do
    end do

    !Next we generate the overal form of the hamiltonian, by adding and normalizing
    allocate(H(num,num))
    H=(kinetic+pot)/Hmax

    !The next step is to begin expanding in time or using the cheby polynomials to do this for us
    allocate(psit(num,num))
    allocate(psishebyexp(num,npoly))
    do j=0,num-1
        psit(j,0)=wave(j)
        psishebyexp(j,0)=wave(j)
    end do

    !Now that we have set up our matrices, we now undergo the recurrence realation for the ChebyShev Polynomials
    allocate(interm(num))
    interm=DOT_PRODUCT(H(:,:),psishebyexp(:,1))
    do j=0,num-1
        psishebyexp(j,1)=interm(j)
    end do
    do n=1, num
        psishebyexp(j,0)=psit(j,n)
        do jj=3,npoly
            interm=2*DOT_PRODUCT(H(:,:),psishebyexp(:,jj-1))
            do j=0,num-1
                psishebyexp(j,jj)=interm(j)
                psishebyexp(j,jj)=psishebyexp(j,jj)-psishebyexp(j,j-2)
            end do
        end do
        do j=3,npoly
            interm=DOT_PRODUCT(psishebyexp(:,:),jntau(:))
            psit(j,n+1)=interm(j)
        end do
    end do
    print*,psit

end program main

    
