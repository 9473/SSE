``` fortran
! subroutine loopupdate()
!     use Configuration
!     implicit none

!     !parameters:
!     integer :: i,n,l,b,op,s1,s2,v0,v1,v2,top,a,g
!     real(8),external :: ran
!     integer,external :: ir
!     integer :: stack(0:4*mm-1)

!     frstspinop = -1
!     lastspinop = -1

!     vertexlist(:) = -2
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

!                 lastspinop(s1)=v0+2
!                 lastspinop(s2)=v0+3
!             elseif (mod(op,4)==3) then !---diag site
!                 g = (op+1)/4

!                 s1 = bsites(1,g)
!                 v1=lastspinop(s1) !use imformation of lastspin(i)


!                 if(v1/=-1) then
!                     vertexlist(v1)=v0
!                     vertexlist(v0)=v1
!                 else 
!                     frstspinop(s1)=v0 !------ leg0
!                 endif

!                 vertexlist(v0+1)=-2
!                 vertexlist(v0+3)=-2

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
!             vertexlist(v2)=v1
!             vertexlist(v1)=v2
!         endif
!     enddo

!     !----------flip loop----------------!
!     do v0=0,4*mm-1,2 !*必须每个leg遍历那样来找       
        
!         if((vertexlist(v0)<0)) cycle !空算符or被访问过
        
!         stack(:) = 0
!         top = 0
!         stack(top)=v0
!         if (ran()<0.5d0) then

!             do while(top>=0)
!                 v1 = stack(top)
!                 top = top-1

!                 if (vertexlist(v1)<0) cycle
!                 v2 = vertexlist(v1)
!                 if ((v2>=0).and.(vertexlist(v2)>=0)) then
!                     top = top + 1
!                     stack(top) = v2
!                 endif
                
!                 vertexlist(v1)=-1
!                 i = opstring(v1/4)
                
!                 if (mod(i,4)==3 .or. mod(i,4)==2) then
!                     a=ieor(mod(i,4),1)!
!                     opstring(v1/4) = opstring(v1/4)-mod(opstring(v1/4),4)+a
!                 endif

!                 if (mod(i,4)==0 .and. i/=0) then
!                     if(vertexlist(ir(v1))>=0) then
!                     !'处理直进'
!                     top = top + 1
!                     stack(top) = ir(v1)
!                     endif
!                     if(vertexlist(ieor(v1,1))>=0) then
!                     !'处理横腿：'
!                     top = top + 1
!                     stack(top) = ieor(v1,1)
!                     endif
!                     if(vertexlist(ieor(ir(v1),1))>=0) then
!                     !'处理斜穿'
!                     top = top + 1
!                     stack(top) = ieor(ir(v1),1)
!                     endif
                    
!                 endif
            
!             enddo
            
!         else
!             do while(top>=0)
!                 v1 = stack(top)
!                 top = top-1

!                 if (vertexlist(v1)<0) cycle 
!                 v2 = vertexlist(v1)
!                 if ((v2>=0).and.(vertexlist(v2)>=0)) then 
!                     top = top + 1
!                     stack(top) = v2
!                 endif
                
!                 vertexlist(v1)=-2
!                 i = opstring(v1/4)

!                 if (mod(i,4)==0 .and. i/=0) then
!                     if(vertexlist(ir(v1))>=0) then
!                     !'处理直进'
!                     top = top + 1
!                     stack(top) = ir(v1)
!                     endif
!                     if(vertexlist(ieor(v1,1))>=0) then
!                     !'处理横腿：'
!                     top = top + 1
!                     stack(top) = ieor(v1,1)
!                     endif
!                     if(vertexlist(ieor(ir(v1),1))>=0) then
!                     !'处理斜穿'
!                     top = top + 1
!                     stack(top) = ieor(ir(v1),1)
!                     endif
                    
!                 endif
!             enddo
!         endif
!     enddo

!     do i=1,nn
!         if (frstspinop(i)/=-1) then
!             if (vertexlist(frstspinop(i))==-1) spin(i)=-spin(i)
!         else
!             if (ran()<0.5) spin(i)=-spin(i)
!         endif
!     enddo

! end subroutine
! ------------------------------!
```
