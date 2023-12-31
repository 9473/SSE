Module Configuration
    ! implicit none
    save

    integer :: lx
    integer :: ly     !----2D Lattice parameter
    integer :: nn     !----2D Lattice numbers of atoms
    integer :: nb     !----numbers of bond
    integer :: nh     !----non-I operators
    integer :: mm     !----cut off

    real(8) :: beta
    real(8) :: aprob
    real(8) :: dprob
    real(8) :: Temperature

    integer, allocatable :: spin(:) !-----restore spin
    integer, allocatable :: bsites(:,:) !---restore No. bonds and sites
    integer, allocatable :: opstring(:) !---restore operators in tau

    integer, allocatable :: frstspinop(:)
    integer, allocatable :: lastspinop(:)
    integer, allocatable :: vertexlist(:)

end Module Configuration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module measurementdata
    ! implicit none
    save

    real(8) :: amag1=0.d0 
    real(8) :: amag2=0.d0 
    real(8) :: data1(2)=0.d0
    real(8) :: data2(2)=0.d0

end module measurementdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





program basic_heisenberg_sse
    use Configuration
    implicit none

    ! parameters:
    integer :: i,j,nbins,msteps,isteps
    


    open(10,file='input.in',status='old')
    read(10,*)lx,ly,Temperature,nbins,msteps,isteps
    close(10)

    beta = 1.0/Temperature
    open(16,file='T.dat',position='append')
    write(16,'((e16.8,2x))') dble(Temperature)
    close(16)


    ! initial all:
    call initran(1)     !----input a seed to get random number!
    call makelattice()  !----constructe the 2D lattice and bond number
    call initconfig()   !----initial configuration

    !use aprob and bprob in diagupdate
    aprob = 0.5d0*beta*nb
    dprob = 1.d0/(0.5d0*beta*nb)

    ! program main
    ! equilibration:
    do i = 1,isteps     !----equilibration need steps
        call diagonalupdate()
        call loopupdate()
        call adjustcutoff(i)
    enddo

    open(10,file='results.txt',status='replace')
    write(10,*)'Finished equilibration, M = ',mm
    close(10)

    ! sample:
    do j = 1,nbins
        do i = 1,msteps
            call diagonalupdate()        
            call loopupdate()
            call measureobservables()
        enddo
        call writeresults(msteps,j)
    enddo


    call deallocateall()


end program


subroutine diagonalupdate()
    use configuration
    implicit none

    !parameters:
    integer :: i,b,op
    real(8) :: p
    real(8),external :: ran

    !-------check bond type
    !----if I, none do; if diag, none do; if off-diag, exchange
    !----if I, change to diag, if diag, change to I
    do i = 0,mm-1
        op = opstring(i)
        if (op==0) then !----I operator
            b = min( int(ran()*nb)+1,nb ) !-----here insert a diag
            if ( spin(bsites(1,b)) /= spin(bsites(2,b)) ) then !---only anti-paralleled can insert!
                ! ran < aprob/(M-n),aprob/(M-n)>=1
                if(aprob>=dfloat(mm-nh) .or. aprob>=ran()*(mm-nh)) then
                    opstring(i) = 2*b
                    nh = nh+1
                endif
            endif
        elseif(mod(op,2)==0) then !----diag
            !ran < dprob*(M-n+1), dprob*(M-n+1)>=1
            p = dprob*(mm-nh+1)
            if(p>=1.d0 .or. p>=ran()) then
                opstring(i) = 0
                nh=nh-1
            endif
        else
            b=op/2
            spin(bsites(1,b))=-spin(bsites(1,b))
            spin(bsites(2,b))=-spin(bsites(2,b))
        endif
    enddo
end subroutine

!--------check type of bond-----------------!
! opstring(i) = op
! opstring = 2b+a-1
! index a = 1, diag; index a = 2, off-diag
! a=1, op=even, diag; a=2, op=odd, off-diag
! b=op/2, op=odd => op = 2b => b=op/2
! b=op/2, op=even => op=2b => b=op/2
!------------------------------!

!--------bond and site--------------------!
! bond(b) left = bond(b) = i
! bond(b) right = bond(b)+1 = j
! bond(b) up = bond(b) + lx = j
! bond(b) down = bond(b) = i
!--------------------------------!

!------exchange I <=> diag--------------!
!----insert diag
! p = beta/2 * Nb / (M-n)
!----delete diag
! p = 2(M-n+1)/Nb/beta

! aprob = nb*beta/2, dprob = 1/(1/2*beta*nb) restore in program main
! ran < aprob/(M-n), ran < dprob*(M-n+1)
!--------------------------------!

!------check spin(i),spin(j)-------------!
! (+1 -1)H(1,b)(+1 -1) = 1/2
! (+1 -1)H(2,b)(-1 +1) = 1/2
!-------------------------------!


subroutine loopupdate()
    use Configuration
    implicit none

    !parameters:
    integer :: i,n,l,b,op,s1,s2,v0,v1,v2
    real(8),external :: ran

    frstspinop = -1
    lastspinop = -1

    !v_0 -> each p
    do v0=0,4*mm-1,4
        op=opstring(v0/4) !--------get op = v_0/4
        if( op/=0 ) then !------if non-I
            b = op/2 !----get bond location
            s1=bsites(1,b) !---get i
            s2=bsites(2,b) !---get j

            v1=lastspinop(s1) !use imformation of lastspin(i)
            v2=lastspinop(s2) !lastspin(j)

            !----should check if i,j has visited-----!
            ! if hasn't checked, v1 must is leg0, v2 = leg1
            ! if has checked, v1 must is old vertex's leg2, v2=leg3
            !-----------------------------------------!

            if(v1/=-1) then
                vertexlist(v1)=v0
                vertexlist(v0)=v1
            else 
                frstspinop(s1)=v0 !------ leg0
            endif

            if(v2/=-1) then
                vertexlist(v2)=v0+1
                vertexlist(v0+1)=v2
            else
                frstspinop(s2)=v0+1 !----- leg1
            endif


            !write now vertex leg2,leg3 in lastspin
            lastspinop(s1)=v0+2
            lastspinop(s2)=v0+3

        else !----------I
            vertexlist(v0:v0+3)=0 !----no vertex
        endif
    enddo

    !periodic boundary condition

    do s1=1,nn
        v1=frstspinop(s1)
        if (v1/=-1) then
            v2=lastspinop(s1)
            vertexlist(v2)=v1
            vertexlist(v1)=v2
        endif
    enddo


    !----------flip loop----------------!
    do v0=0,4*mm-1,2       
        if (vertexlist(v0)<1) cycle
        v1=v0
        if (ran()<0.5d0) then
            do 
                opstring(v1/4)=ieor(opstring(v1/4),1)
                vertexlist(v1)=-1
                v2=ieor(v1,1)
                v1=vertexlist(v2)
                vertexlist(v2)=-1
                if (v1==v0) exit
            enddo
        else
            do 
                vertexlist(v1)=0
                v2=ieor(v1,1)
                v1=vertexlist(v2)
                vertexlist(v2)=0
                if (v1==v0) exit
            enddo
        endif
    enddo

    do i=1,nn
        if (frstspinop(i)/=-1) then
            if (vertexlist(frstspinop(i))==-1) spin(i)=-spin(i)
        else
            if (ran()<0.5) spin(i)=-spin(i)
        endif
    enddo

end subroutine

!------vertex------------!
! v(p) = 4p+i, i =0,1,2,3
! v_0(p) = 4p ---> op = v_0/4
! each v0 has only two cases:
! --1. meet leg0
! --2. meet leg2
!------------------------!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  measure 

subroutine measureobservables()
    use Configuration
    use measurementdata
    implicit none

    integer :: i,b,op,s1,s2,am
    real(8) :: am1,am2
   
    am=0
    do i=1,nn
       am=am+spin(i)*(-1)**(mod(i-1,lx)+(i-1)/lx) !staggered mag
    enddo      
    am=am/2 !m_s = 1/2 sum_i phi_i sigma_i
    am1=0.d0
    am2=0.d0
    do i=0,mm-1
       op=opstring(i)
       if (op==0) then
           cycle
       elseif (mod(op,2)==1) then        
          b=op/2
          s1=bsites(1,b)
          s2=bsites(2,b)
          spin(s1)=-spin(s1)
          spin(s2)=-spin(s2)
          am=am+2*spin(s1)*(-1)**(mod(s1-1,lx)+(s1-1)/lx) !两个自旋调换,差值是两倍
       endif
       am1=am1+dfloat(abs(am)) !交错磁化的绝对值 |m|
       am2=am2+dfloat(am)**2 !交错磁化的平方 m^2
    enddo
    if (nh/=0) then
       am1=am1/nh
       am2=am2/nh
    else
       am1=dfloat(abs(am))
       am2=dfloat(am)**2
    endif

    amag1=amag1+am1 !每个ns下的值传递给amag
    amag2=amag2+am2

end subroutine


subroutine writeresults(msteps,bins)
    use Configuration
    use measurementdata
    implicit none
    
    integer :: i,msteps,bins
    real(8) :: wdata1(2),wdata2(2)
   
    amag1=amag1/msteps
    amag2=amag2/msteps

    amag1=amag1/nn !计算平均磁化，所以要除以格点数
    amag2=amag2/nn/nn

    data1(1)=data1(1)+amag1 !每个bin下的值，累积传递给data1
    data1(2)=data1(2)+amag2

   
    data2(1)=data2(1)+amag1**2 !物理量的平方传递给data2
    data2(2)=data2(2)+amag2**2

    do i=1,2
        wdata1(i)=data1(i)/bins !对bin求平均值传递给wdata
        wdata2(i)=data2(i)/bins
        wdata2(i)=sqrt(abs(wdata2(i)-wdata1(i)**2)/bins) !计算物理量的方差:sqrt[ (<A^2> - <A>^2) /bins]
     enddo

    ! open(10,file='results.txt',status='replace')
    ! write(10,*)' Cut-off L : ',mm
    ! write(10,*)' Number of bins completed : ',bins

    ! write(10,*)' ========================================='
    ! write(10,10)' -E/N       : ',wdata1
    ! write(10,*)' ========================================='
    ! 10 format(1x,a,2f14.8)
    ! close(10)


    ! open(16,file='T.dat',position='append')
    ! write(16,'((e16.8,2x))') dble(Temperature)
    ! close(16)

    open(16,file='absm.dat',position='append')
    write(16,'((e16.8,2x))') dble(wdata1(1))
    close(16)

    open(15,file='m2.dat',position='append')
    write(15,'((e16.8,2x))') dble(wdata1(2))
    close(15)

    ! open(16,file='absmerror.dat',position='append')
    ! write(16,'((e16.8,2x))') dble(wdata2(1))
    ! close(16)

    ! open(16,file='m2error.dat',position='append')
    ! write(16,'((e16.8,2x))') dble(wdata2(2))
    ! close(16)



    amag1=0.d0
    amag2=0.d0


end subroutine

! measure ends
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine adjustcutoff(step)

    use configuration
    implicit none

    integer, allocatable :: stringcopy(:)
    integer :: mmnew,step

    mmnew=nh+nh/3
    if (mmnew<=mm) return !M is enough and doesn't change

    allocate(stringcopy(0:mm-1))
    stringcopy(:)=opstring(:)
    deallocate(opstring)
    allocate(opstring(0:mmnew-1))
    opstring(0:mm-1)=stringcopy(:)
    opstring(mm:mmnew-1)=0
    deallocate(stringcopy)

    mm=mmnew
    deallocate (vertexlist)
    allocate(vertexlist(0:4*mm-1))

    open(unit=10,file='results.txt',status='replace')
    write(10,*)' Step: ',step,'  Cut-off L: ',mm
    close(10)
    
end subroutine adjustcutoff

subroutine initconfig()
    use configuration
    implicit none

    integer :: i
    real(8),external :: ran

    allocate(spin(nn))
    do i=1,nn
    spin(i)=2*int(2.*ran())-1
    enddo

    mm=20
    allocate(opstring(0:mm-1))
    opstring(:)=0
    nh=0

    allocate(frstspinop(nn))
    allocate(lastspinop(nn))
    allocate(vertexlist(0:4*mm-1))
    
end subroutine initconfig

subroutine makelattice()
    use configuration
    implicit none
    integer :: s,x1,y1,x2,y2

    nn = lx*ly
    nb = 2*nn
    allocate(bsites(2,nb))

    do y1 = 0,ly-1 !! start from 0!
    do x1 = 0,lx-1
        !number:
        s = y1*lx + x1 + 1

        ! x---->
        x2 = mod(x1+1,lx) !x+1
        y2 = y1
        
        bsites(1,s) = s   !bond left = bond(b)
        bsites(2,s) = 1+x2+y2*lx !bond right = bond(b)+1

        ! y---->
        x2 = x1
        y2 = mod(y1+1,ly) !y+1

        bsites(1,s+nn) = s !bond down = bond(b)
        bsites(2,s+nn) = 1+x2+y2*lx !bond up = bond(b)+lx
    enddo
    enddo

end subroutine makelattice
!--------bond and site--------------------!
! bond(b) left = bond(b) = i
! bond(b) right = bond(b)+1 = j
! bond(b) up = bond(b) + lx = j
! bond(b) down = bond(b) = i
!--------------------------------!


subroutine deallocateall()
    use configuration
    implicit none

    deallocate (spin)
    deallocate (bsites)
    deallocate (opstring)
    deallocate (frstspinop)
    deallocate (lastspinop)
    deallocate (vertexlist)
    
end subroutine deallocateall


real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!
