# J1J2 model for SSE
和标准的Heisenberg model没有太大的差别，如果我们同样假设加入的常数使得对角bond算符插入在自旋相同的两个点上的权重为零，这样的话对角bond算符只能插入到自旋反平行的两个格点之间。不过注意到，在奇数还是偶数上插入的是不同的对角算符.  

$$
H_{diag,J_1} = \frac{J_1}{4}-J_1S^zS^z
$$

$$
H_{diag,J_2} = \frac{J_2}{4}-J_2S^zS^z
$$

而且这种方式会使得对角bond翻转到各自的非对角bond的时候权重不变，因此可以loop翻转前后权重不变.  

------

以下更新方式符合在x轴上奇数键插入J1，在偶数键上插入J2，其余键（竖着的键）均插入J2.
<img width="253" alt="image" src="https://github.com/9473/SSE/assets/59651278/c37c0fb1-0bf3-4a83-86d5-40c62f0a2e8c">

此计数方式更详细见SSE教程. 可以看到，横键的范围为 1至nn, 竖键的范围为 nn+1至nb(其中nb=2nn)

### 对角更新

Set $J_1 = 1$,  $jr = \frac{J_2}{J_1}$. 那么对角更新的概率:

$p_{ist} = \frac{\beta N_b \cdot 1/2}{M-n}$.  可以认为是插入 $J_1$ 对角bond的概率，那么插入 $J_2$ 对角bond的概率可以很简单地认为： $p_{ist} = \frac{\beta N_b \cdot 1/2}{M-n}*J_r$.  在这样的情况下，插入 $J_1$ 对角bond的概率可以认为是 $p_{ist J2}/jr$.

 相应的，删去算符在 $J_2$ bond的情况下也要做相应的 $p*J_2$ 的修改。

```fortran
! in main:
aprob=0.5d0*beta*nb
dprob=1.d0/(0.5d0*beta*nb)
!提前定义好全局变量 jr

!-----------------!
subroutine diagonalupdate()
 use configuration; implicit none

 integer :: i,b,op,s
 real(8) :: p
 real(8), external :: ran

 do i=0,mm-1
    op=opstring(i)
    if (op==0) then       
       b=min(int(ran()*nb)+1,nb)
       if (spin(bsites(1,b))/=spin(bsites(2,b))) then
	  p=aprob  !默认插入J2
	    if (b<=nn.and.mod((b-1),2)==0) p=p/jr !表明横轴奇数键上插入J1
          if (p>=dfloat(mm-nh).or.p>=ran()*(mm-nh)) then
             opstring(i)=2*b
             nh=nh+1 
          endif
       endif
    elseif (mod(op,2)==0) then  !说明这里对角算符的计数方式是2的偶数倍
       b=op/2        
       p=dprob*(mm-nh+1)
		if (b<=nn.and.mod((b-1),2)==0) p=p*jr
       if (p>=1.d0.or.p>=ran()) then
          opstring(i)=0
          nh=nh-1
       endif
    else
       b=op/2
       spin(bsites(1,b))=-spin(bsites(1,b))
       spin(bsites(2,b))=-spin(bsites(2,b))
    endif
 enddo

end subroutine diagonalupdate
```

对角更新修改完成.

------
### Loop更新

和原来的一样.  
```fortran
subroutine loopupdate()
!-----------------------!
 use configuration; implicit none

 integer :: i,n,l,b,op,s1,s2,v0,v1,v2,s
 real(8), external :: ran

 frstspinop(:)=-1
 lastspinop(:)=-1

 do v0=0,4*mm-1,4
    op=opstring(v0/4)
    if (op/=0) then
       b=op/2
       s1=bsites(1,b)
       s2=bsites(2,b)
       v1=lastspinop(s1)
       v2=lastspinop(s2)
       if (v1/=-1) then
          vertexlist(v1)=v0
          vertexlist(v0)=v1
       else
          frstspinop(s1)=v0
       endif
       if (v2/=-1) then
          vertexlist(v2)=v0+1
          vertexlist(v0+1)=v2
       else
          frstspinop(s2)=v0+1
       endif
       lastspinop(s1)=v0+2
       lastspinop(s2)=v0+3
    else
       vertexlist(v0:v0+3)=0
    endif
 enddo
 do s1=1,nn
    v1=frstspinop(s1)
    if (v1/=-1) then
        v2=lastspinop(s1)
        vertexlist(v2)=v1
        vertexlist(v1)=v2
    endif
 enddo

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

end subroutine loopupdate
```

Loop更新结束

------
### make lattice
J1J2 model(若)仍然使用方晶格格子，则代码一致，不会因为J1,J2参数造成任何差异.  

```fortran
!这里演示的是PBC 2D lattice
subroutine makelattice()
!------------------------!
 use configuration; implicit none

 integer :: s,x1,x2,y1,y2,i,x,y

 nn=lx*ly
 nb=2*nn 

 allocate(bsites(2,nb))

 do y1=0,ly-1
 do x1=0,lx-1
    s=1+x1+y1*lx
    x2=mod(x1+1,lx)
    y2=y1
    bsites(1,s)=s
    bsites(2,s)=1+x2+y2*lx
    x2=x1
    y2=mod(y1+1,ly)
    bsites(1,s+nn)=s
    bsites(2,s+nn)=1+x2+y2*lx       
 enddo
 enddo
end subroutine makelattice


------

另外一种处理J1J2 model的方式是在造格子上对所属J1和J2的格点进行改造，但在这里就不展示了。
