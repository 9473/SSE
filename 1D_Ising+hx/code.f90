Module Configuration
    ! implicit none
    save

    integer :: lx
    integer :: ly     !----2D Lattice parameter
    integer :: nn     !----2D Lattice numbers of atoms
    integer :: nb     !----numbers of bond
    integer :: nh,nbond,nsite,nhh    !----non-I operators
    integer :: mm     !----cut off

    real(8) :: beta
    real(8) :: aprob,apb
    real(8) :: dprob,dprob_s,dprob_b
    real(8) :: h
    real(8) :: coupling

    integer, allocatable :: spin(:) !-----restore spin
    integer, allocatable :: bsites(:,:) !---restore No. bonds and sites
    integer, allocatable :: opstring(:) !---restore operators in tau

    integer, allocatable :: frstspinop(:)
    integer, allocatable :: lastspinop(:)
    integer, allocatable :: vertexlist(:)
    ! integer, allocatable :: stack(:)
    ! integer :: top,zz
    ! integer, allocatable :: visited(:)

end Module Configuration


module measurementdata
    ! implicit none
    save

    real(8) :: enrg1=0.d0
    real(8) :: ener =0.d0
    real(8) :: enrg2=0.d0
    real(8) :: amag1=0.d0 
    real(8) :: amag2=0.d0 
    real(8) :: asusc=0.d0 
    real(8) :: stiff=0.d0 
    real(8) :: ususc=0.d0
    real(8) :: data1(8)=0.d0
    real(8) :: data2(8)=0.d0

end module measurementdata



program TFIS_sse
    use Configuration
    implicit none

    ! parameters:
    integer :: i,j,nbins,msteps,isteps
    

    open(10,file='input.in',status='old')
    read(10,*)lx,nbins,msteps,isteps
    close(10)

    beta = 20
    h=0.8
    coupling = 1.0
    ! open(16,file='T.dat',position='append')
    ! write(16,'((e16.8,2x))') dble(Temperature)
    ! close(16)


    ! initial all:
    call initran(1)     !----input a seed to get random number!
    call makelattice()  !----constructe the 2D lattice and bond number
    call initconfig()   !----initial configuration

    !use aprob and bprob in diagupdate
    aprob = (2.0*nb*coupling+nn*h)*beta !(2JN_b + hN)beta
        apb = (2.0*nb*coupling)/(2.0*nb*coupling+nn*h)
    dprob = 1.0/aprob

    
    ! program main
    ! equilibration:
    do i = 1,isteps     !----equilibration need steps
        call diagonalupdate()
        call loopupdate()
        call adjustcutoff(i)
    enddo
    print*,nh,"<----nh",mm,"<----M"
    ! print*,opstring(:)
    ! print*,opstring(2),opstring(4),opstring(6)


    ! open(10,file='results.txt',status='replace')
    ! write(10,*)'Finished equilibration, M = ',mm
    ! close(10)
    

    ! ! sample:
    do j = 1,nbins
        do i = 1,msteps
            call diagonalupdate()        
            call loopupdate()
            call measureobservables()
        enddo
        call writeresults(msteps,j)
    enddo
    ! ! print*,nh,"<----nh",mm,"<----M"
    ! print*,"vertexlist:"
    ! print*,vertexlist(:)

    call deallocateall()

    print*,"main end"

end program

subroutine diagonalupdate()
    use configuration
    implicit none

    !parameters:
    integer :: i,b,op,g
    real(8) :: p
    real(8),external :: ran

  !--!*(1/17)重新写了这部分的程序，主要是更换概率发生先后的逻辑--------------!
  !opstring(i) = op, op = 4*b + a - 1
  !H_{-1,g} = off-diag site: op = 4*g -2, mod(op,4) = 2
  !H_{0,g} = diag site: op = 4*g -1, mod(op,4) = 3
  !H_{1,b} = diga bond: op = 4*b, mod(op,4) = 0,op/=0
  !-------------------------------------------------------------------!
    do i = 0,mm-1
        op = opstring(i)
        if (op==0) then !----I operator
            !now need to insert operator:
            
            if (apb>=ran()) then! insert bond
                if ( aprob>=dfloat(mm-nh) .or. aprob>=ran()*(mm-nh) ) then!(beta/(M-n)) *(2JN_b + hN) 是否要插键的概率
                !sure insert bond
                
                    b = min( int(ran()*nb)+1,nb ) !-----here insert a diag
                    ! print*,"b=",b
                    if (spin(bsites(1,b))==spin(bsites(2,b))) then !只有平行矩阵元起作用 !反平行自旋权重为0
                        opstring(i) = 4*b
                        nh = nh+1
                        ! print*,"insert bond diag"
                    ! else
                    ! print*,"site:",spin(:)
                    endif
                endif
            else !else insert diag site

                if ( aprob>=dfloat(mm-nh) .or. aprob>=ran()*(mm-nh) ) then
                    g = min( int(ran()*nn)+1,nn ) !diagonal site
                    opstring(i) = 4*g - 1 !4b+a-1,a=0
                    nh = nh+1
                    ! print*,"insert site diag"
                endif
            endif !不符合插键概率，不插直接结束

        elseif(((op/=0).and.(mod(op,4)==0)).or.(mod(op,4)==3)) then !----bond-diag or g-diag !* 无论遇到哪种，删去的概率相同，因此写在一起
            !ran < dprob*(M-n+1), dprob*(M-n+1)>=1
            p = dprob*(mm-nh+1)
            if(p>=1.d0 .or. p>=ran()) then
                opstring(i) = 0
                nh=nh-1
                ! print*,"delete a diag"
            endif
        else
            if(mod(op,4)==2) then
                g=(op+2)/4
                spin(bsites(1,g))=-spin(bsites(1,g))
            endif
        endif
    enddo

end subroutine


subroutine loopupdate()
    use Configuration
    implicit none

    !parameters:
    integer :: i,n,l,b,op,s1,s2,v0,v1,v2,top,a,g
    real(8),external :: ran
    integer,external :: ir
    integer :: stack(0:8*mm-1)

    frstspinop = -1
    lastspinop = -1

    vertexlist(:) = -2
    !v_0 -> each p
    do v0=0,4*mm-1,4
        op=opstring(v0/4) !--------get op = opstring(i), i = p = v0/4
        if( op/=0 ) then !------if non-I
            ! print*,"op=",op
            if (mod(op,4)==0) then !---diag bond
            ! print*,"it's diag bond,     ","mod(op,4)=",mod(op,4)
                b = op/4 !----get bond location
                
                s1=bsites(1,b) !---get i
                s2=bsites(2,b) !---get j
                v1=lastspinop(s1) !use imformation of lastspin(i)
                v2=lastspinop(s2) !lastspin(j)

                if(v1/=-1) then
                    ! print*,v1,"<---v1",v0,"<---v0"
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

                lastspinop(s1)=v0+2
                lastspinop(s2)=v0+3
            elseif (mod(op,4)==3) then !---diag site
                g = (op+1)/4

                s1 = bsites(1,g)
                v1=lastspinop(s1) !use imformation of lastspin(i)


                if(v1/=-1) then
                    vertexlist(v1)=v0
                    vertexlist(v0)=v1
                else 
                    frstspinop(s1)=v0 !------ leg0
                endif

                vertexlist(v0+1)=-2
                vertexlist(v0+3)=-2

                lastspinop(s1)=v0+2

            else!*非对角算符:
                if (mod(op,4)==2) then !---off-diag site

                    g = (op+2)/4
                    s1 = bsites(1,g)
                    v1=lastspinop(s1) !use imformation of lastspin(i)

                if(v1/=-1) then
                    ! print*,v1,"<---v1",v0,"<---v0"
                    vertexlist(v1)=v0
                    vertexlist(v0)=v1
                else 
                    frstspinop(s1)=v0 !------ leg0
                endif
                    vertexlist(v0+1)=-2
                    vertexlist(v0+3)=-2
                    lastspinop(s1)=v0+2
                
                endif !for off-diag site
            endif !for check mod(op,4)
        else !----------I
            vertexlist(v0:v0+3)=-2 !----no vertex
        endif
    enddo
    

    !periodic boundary condition
    do s1=1,nn
        v1=frstspinop(s1)
        if ((v1/=-1)) then
            v2=lastspinop(s1)
            vertexlist(v2)=v1
            vertexlist(v1)=v2
        endif
    enddo

    !----------flip loop----------------!
    do v0=0,4*mm-1,2 !*必须每个leg遍历那样来找       
        
        if((vertexlist(v0)<0)) cycle !空算符or被访问过
        
        stack(:) = 0
        top = 0
        stack(top)=v0
        if (ran()<0.5d0) then

            do while(top>=0)
                v1 = stack(top)
                top = top-1

                if (vertexlist(v1)<0) cycle
                v2 = vertexlist(v1)
                if ((v2>=0).and.(vertexlist(v2)>=0)) then
                    top = top + 1
                    stack(top) = v2
                endif
                
                vertexlist(v1)=-1
                i = opstring(v1/4)
                
                if (mod(i,4)==3 .or. mod(i,4)==2) then
                    a=ieor(mod(i,4),1)!
                    opstring(v1/4) = opstring(v1/4)-mod(opstring(v1/4),4)+a
                endif

                if (mod(i,4)==0 .and. i/=0) then
                    if(vertexlist(ir(v1))>=0) then
                    !'处理直进'
                    top = top + 1
                    stack(top) = ir(v1)
                    endif
                    if(vertexlist(ieor(v1,1))>=0) then
                    !'处理横腿：'
                    top = top + 1
                    stack(top) = ieor(v1,1)
                    endif
                    if(vertexlist(ieor(ir(v1),1))>=0) then
                    !'处理斜穿'
                    top = top + 1
                    stack(top) = ieor(ir(v1),1)
                    endif
                    
                endif
            
            enddo
            
        else
            do while(top>=0)
                v1 = stack(top)
                top = top-1

                if (vertexlist(v1)<0) cycle 
                v2 = vertexlist(v1)
                if ((v2>=0).and.(vertexlist(v2)>=0)) then 
                    top = top + 1
                    stack(top) = v2
                endif
                
                vertexlist(v1)=-2
                i = opstring(v1/4)

                if (mod(i,4)==0 .and. i/=0) then
                    if(vertexlist(ir(v1))>=0) then
                    !'处理直进'
                    top = top + 1
                    stack(top) = ir(v1)
                    endif
                    if(vertexlist(ieor(v1,1))>=0) then
                    !'处理横腿：'
                    top = top + 1
                    stack(top) = ieor(v1,1)
                    endif
                    if(vertexlist(ieor(ir(v1),1))>=0) then
                    !'处理斜穿'
                    top = top + 1
                    stack(top) = ieor(ir(v1),1)
                    endif
                    
                endif
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
! ------------------------------!


subroutine measureobservables()
    use Configuration
    use measurementdata
    implicit none

    integer :: i,b,op,s1,s2,am,jj(0:1)
    real(8) :: am1,am2,ax1
   
    ! am=0
    ! do i=1,nn
    !    am=am+spin(i)*(-1)**(mod(i-1,lx)+(i-1)/lx)
    ! enddo      
    ! am=am/2
    ! am1=0.d0
    ! am2=0.d0
    ! ax1=0.d0
    ! jj(:)=0
    ! do i=0,mm-1
    !    op=opstring(i)
    !    if (op==0) then
    !        cycle
    !    elseif (mod(op,2)==1) then        
    !       b=op/2
    !       s1=bsites(1,b)
    !       s2=bsites(2,b)
    !       spin(s1)=-spin(s1)
    !       spin(s2)=-spin(s2)
    !       jj((b-1)/nn)=jj((b-1)/nn)+spin(s2)
    !       am=am+2*spin(s1)*(-1)**(mod(s1-1,lx)+(s1-1)/lx)
    !    endif
    !    ax1=ax1+dfloat(am)
    !    am1=am1+dfloat(abs(am))
    !    am2=am2+dfloat(am)**2
    ! enddo
    ! if (nh/=0) then
    !    ax1=(ax1**2+am2)/(dfloat(nh)*dfloat(nh+1))
    !    am1=am1/nh
    !    am2=am2/nh
    ! else
    !    am1=dfloat(abs(am))
    !    am2=dfloat(am)**2
    !    ax1=am2
    ! endif
   

    nsite=0
    do i = 0,mm-1,1
        if(mod(opstring(i),4)==3 .or. mod(opstring(i),4)==2) nsite =nsite+1
    enddo
    nhh=0
    do i = 0,mm-1,1
        if(mod(opstring(i),4)==3) nhh =nhh+1
    enddo
    
    enrg1=enrg1+dfloat(nh) 
    ener =ener +dfloat(nh-nsite) !This is not the right energy
    enrg2=enrg2+dfloat(nh-nhh) 
    amag1=amag1+am1
    amag2=amag2+am2
    asusc=asusc+ax1
    stiff=stiff+0.5d0*(dfloat(jj(0))**2+dfloat(jj(1))**2)
    ususc=ususc+dfloat(sum(spin)/2)**2
   

end subroutine


subroutine writeresults(msteps,bins)
    use Configuration
    use measurementdata
    implicit none
    
    integer :: i,msteps,bins
    real(8) :: wdata1(8),wdata2(8)
   
    enrg1=enrg1/msteps
    ener =ener /msteps
    enrg2=enrg2/msteps
    amag1=amag1/msteps
    amag2=amag2/msteps
    asusc=asusc/msteps
    stiff=stiff/msteps
    ususc=ususc/msteps
   
    enrg2=enrg2 /(-beta)/nn + coupling 
    enrg1=enrg1/(beta)/nn-coupling-h
    ener =ener /(-beta)/nn + coupling! + h!*对常数项进行处理
    amag1=amag1/nn
    amag2=amag2/nn
    asusc=beta*asusc/nn    
    ususc=beta*ususc/nn
    stiff=stiff/(beta*nn)
   
    data1(1)=data1(1)+enrg1
    data1(2)=data1(2)+enrg2
    data1(3)=data1(3)+amag1
    data1(4)=data1(4)+amag2
    data1(5)=data1(5)+asusc
    data1(6)=data1(6)+stiff
    data1(7)=data1(7)+ususc
    data1(8)=data1(8)+ener
   
    data2(1)=data2(1)+enrg1**2
    data2(2)=data2(2)+enrg2**2
    data2(3)=data2(3)+amag1**2
    data2(4)=data2(4)+amag2**2
    data2(5)=data2(5)+asusc**2
    data2(6)=data2(6)+stiff**2
    data2(7)=data2(7)+ususc**2
    data2(8)=data2(8)+ener **2
   
    do i=1,8
       wdata1(i)=data1(i)/bins
       wdata2(i)=data2(i)/bins
       wdata2(i)=sqrt(abs(wdata2(i)-wdata1(i)**2)/bins)
    enddo
    
    write(*,*) -wdata1(1),wdata1(2),wdata1(8)
    print*,mm,nh,nsite,nhh
    ! open(10,file='results.txt',status='replace')
    ! write(10,*)' Cut-off L : ',mm
    ! write(10,*)' Number of bins completed : ',bins
    ! write(10,*)' ========================================='
    ! write(10,10)' -E/N       : ',wdata1(1),wdata2(1)
    ! write(10,10)'  C/N       : ',wdata1(2),wdata2(2)
    ! write(10,10)'  <|m|>     : ',wdata1(3),wdata2(3)
    ! write(10,10)'  S(pi,pi)  : ',wdata1(4),wdata2(4)
    ! write(10,10)'  X(pi,pi)  : ',wdata1(5),wdata2(5)
    ! write(10,10)'  rho_s     : ',wdata1(6),wdata2(6)
    ! write(10,10)'  X(0,0)    : ',wdata1(7),wdata2(7)
    ! write(10,*)' ========================================='
    ! 10 format(1x,a,2f14.8)
    ! close(10)
   
    enrg1=0.d0
    enrg2=0.d0
    amag1=0.d0
    amag2=0.d0
    asusc=0.d0
    stiff=0.d0
    ususc=0.d0
    ener =0.d0

end subroutine


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

    mm=20 !20
    allocate(opstring(0:mm-1))
    opstring(:)=0
    nh=0
    nhh=0

    allocate(frstspinop(nn))
    allocate(lastspinop(nn))
    allocate(vertexlist(0:4*mm-1))
    ! allocate(stack(0:8*mm-1))

    
end subroutine initconfig

subroutine makelattice()
    use configuration
    implicit none
    integer :: s,x1,y1,x2,y2

    nn = lx
    nb = nn
    allocate(bsites(2,nb))

    do x1 = 0,lx-1
        !number:
        s = x1 + 1
        ! x---->
        x2 = mod(x1+1,lx) !x+1
        
        bsites(1,s) = s   !bond left = bond(b)
        bsites(2,s) = 1+x2

    enddo

end subroutine makelattice


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


integer function ir(x)
 implicit none
integer :: x,x1,x2
x1=x/4
x2=x-x1*4
ir=mod(x2+2,4)+4*x1
return
end function ir !这个子程序确保能够找到vertex x的对面点


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


!----以下是另一种更新方式，在这里备份，需要用的时候
!---将stack设为全局变量，并同时在初始化、adjustcuttoff和deallocate中设置stack
!----将top，zz设为全局变量
!  subroutine copyloopupdate()
! !-----------------------!
!     use Configuration
!     implicit none
!     !parameters:
!     integer :: i,n,l,b,op,s1,s2,v0,v1,v2,a,g
!     real(8),external :: ran
!     integer,external :: ir

!     frstspinop = -1
!     lastspinop = -1

!     vertexlist(:) = 0
!     !v_0 -> each p
!     do v0=0,4*mm-1,4
!         op=opstring(v0/4) !--------get op = opstring(i), i = p = v0/4
!         if( op/=0 ) then !------if non-I
!             ! print*,"op=",op
!             if (mod(op,4)==0) then !---diag bond
!             ! print*,"it's diag bond,     ","mod(op,4)=",mod(op,4)
!                 b = op/4 !----get bond location
                
!                 s1=bsites(1,b) !---get i
!                 s2=bsites(2,b) !---get j
!                 v1=lastspinop(s1) !use imformation of lastspin(i)
!                 v2=lastspinop(s2) !lastspin(j)

!                 if(v1/=-1) then
!                     ! print*,v1,"<---v1",v0,"<---v0"
!                     vertexlist(v1)=v0
!                     vertexlist(v0)=v1
!                 else 
!                     frstspinop(s1)=v0 !------ leg0
!                 endif

!                 if(v2/=-1) then
!                     vertexlist(v2)=v0+1
!                     vertexlist(v0+1)=v2
!                 else
!                     frstspinop(s2)=v0+1 !----- leg1
!                 endif
!                 !write now vertex leg2,leg3 in lastspin
!                 lastspinop(s1)=v0+2
!                 lastspinop(s2)=v0+3
!             elseif (mod(op,4)==3) then !---diag site
!                 g = (op+1)/4
!                 ! s1 = g
!                 ! s2 = g+1!-------!周期性条件在这里是对的吗？——可能还真有问题

!                 s1 = bsites(1,g)
!                 v1=lastspinop(s1) !use imformation of lastspin(i)
!                 ! v2=lastspinop(s2) !lastspin(j)!*可以不从上一个算符那里提取s2

!                 if(v1/=-1) then
!                 ! print*,v1,"<---v1",v0,"<---v0"
!                     vertexlist(v1)=v0
!                     vertexlist(v0)=v1
!                 else 
!                     frstspinop(s1)=v0 !------ leg0
!                 endif
!                     vertexlist(v0+1)=-2
!                     vertexlist(v0+3)=-2
!                 lastspinop(s1)=v0+2
!             else!*非对角算符:
!                 if (mod(op,4)==2) then !---off-diag site

!                     g = (op+2)/4
!                     s1 = bsites(1,g)
!                     v1=lastspinop(s1) !use imformation of lastspin(i)

!                 if(v1/=-1) then
!                     ! print*,v1,"<---v1",v0,"<---v0"
!                     vertexlist(v1)=v0
!                     vertexlist(v0)=v1
!                 else 
!                     frstspinop(s1)=v0 !------ leg0
!                 endif
!                     vertexlist(v0+1)=-2
!                     vertexlist(v0+3)=-2
!                     lastspinop(s1)=v0+2
                
!                 endif !for off-diag site
!             endif !for check mod(op,4)
!         else !----------I
!             vertexlist(v0:v0+3)=-2 !----no vertex
!         endif
!     enddo
    

!     !periodic boundary condition
!     do s1=1,nn
!         v1=frstspinop(s1)
!         if ((v1/=-1)) then
!             v2=lastspinop(s1)
!             ! print*,v1,"<---v1",v2,"<---v2"
!             vertexlist(v2)=v1
!             vertexlist(v1)=v2
!         endif
!     enddo
 

!     top=0
!     do v0=0,4*mm-1,2
!         if (vertexlist(v0)<1) cycle!already in loop
!         ! print*,"where begin"
!         zz=-2
!         if (ran()>0.5) zz=-1
!         call input(v0)
!         !i=0
!         do while(top>=0)
!         ! print*,"where makeloop",top,"<---top"
!             call makeloop()
!             ! if (top==0) exit
!         enddo
!     enddo


!     do i=1,nn
!         if (frstspinop(i)/=-1) then
!             if (vertexlist(frstspinop(i))==-1) spin(i)=-spin(i)
!         else
!             if (ran()<0.5) spin(i)=-spin(i)
!         endif
!     enddo

!  end subroutine

! subroutine input(x)
!  use configuration; implicit none
!  integer :: x

!  top=top+1
!  stack(top)=x

!  end subroutine input

! integer function output()
! use configuration; implicit none

! output=stack(top)
! top=top-1

! end function output


!  subroutine makeloop()
!  use configuration;implicit none
!  integer :: v1,v2,v3,v4,v,ir,i,a
!  integer,external :: output

!  v1=output()
!  v=vertexlist(v1)

!  !!这句是关键，没有这句能量不对!
! !  if ((v>=0).and.(vertexlist(v)/=zz)) call input(v)
!  !!这句是关键，没有这句能量不对!

! !  print*,"where check v and v1" 
!  if (vertexlist(v1)<0) return
!  if ((v>=0).and.(vertexlist(v)/=zz)) call input(v)!*移到了这里
!  vertexlist(v1)=zz
!  i=v1/4

!     if(zz==-1) then !符合概率翻转
!     if  (mod(opstring(i),4)==3 .or. mod(opstring(i),4)==2) then
!         a=ieor(mod(opstring(i),4),1)!
!         opstring(v1/4) = opstring(v1/4)-mod(opstring(v1/4),4)+a
!     endif
!     endif

!     ! print*,"where add bond"
!     if ((mod(opstring(i),4)==0).and.(opstring(i)/=0))  then
!         v2=(ieor(v1,1))
!         v3=(ir(v1))
!         v4=(ieor(v3,1))
!         if (vertexlist(v2)>0) call input(v2)
!         if (vertexlist(v3)>0) call input(v3)
!         if (vertexlist(v4)>0) call input(v4)
!     endif

!  end subroutine makeloop
