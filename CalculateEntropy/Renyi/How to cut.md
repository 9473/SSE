## 纠缠熵的切割方式

这里主要讲的是corners的切割:  

```fortran
alx,aly,lx,ly,cutx,cuty in global

al=alx*aly;

integer, allocatable :: cutsites(:) 
integer, allocatable :: belongA(:) 
integer, allocatable :: inversA(:) 

in initconfig:
allocate(joinorsplit(al))
 joinorsplit(:)=0 This is for B
 
 allocate(cutsites(al))
 allocate(belongA(nn))
 allocate(inversA(nn))

 belongA(:)=0
 inversA(:)=0
 !------------------------------!
 !----cut out part A within-----!
 !------------------------------!
 flag=1
 do j=1,aly
    do i=1,alx
       s=cutx+i
       s1=cuty+j-1
       if(s>lx) then
          s=s-lx
       endif
       if(s1>ly) then
          s1=s1-ly
       endif
       cutsites(flag)=s1*lx+s
       belongA(s1*lx+s)=1       !! belongA =1,in A; =0 in Abar.
       inversA(s1*lx+s)=flag
       flag=flag+1
    enddo
 enddo
```



从 cutx 处横向出发.

cutx+1, cutx+2, .... cutx+alx

从 cuty 纵向出发

Cuty(就在cutx那一排), cuty+1, cuty+2, ... cuty+aly-1.

```bash
s1
s1
s1
s1
s1(s) s s s s
```

第一个点:

```bash
(cutx+1,cuty) 计数方式: cutx+1+(cuty*Lx)
cutsites(1)=cutx+1+(cuty*Lx)
belongA(cutx+1+(cuty*Lx))=1
inversA(cutx+1+(cuty*Lx))=1
```

第二个点：

```bash
(cutx+2,cuty) 计数方式: cutx+2+(cuty*Lx)
cutsites(2)=cutx+2+(cuty*Lx)
belongA(cutx+2+(cuty*Lx))=1
inversA(cutx+2+(cuty*Lx))=2
```

... 

第 alx+1 个点:

```bash
(cutx+1,cuty+1) 计数方式: cutx+1+(cuty+1*Ly)
cutsites(alx+1)=cutx+1+(cuty+1*Ly)
belongA(cutx+1+(cuty+1*Ly))=1
inversA(cutx+1+(cuty+1*Ly))=alx+1
```

显然，这种计数方式与bsites的标记一致.  但是只标记了从cutx+1 to alx, 从cuty to aly-1的部分  

cutsites就像是一个数组的Index用来储存cut那一部分的点的坐标, 其坐标计数方式与bsites一模一样.  而这个坐标必须通过cutsites(index)来获取，这个index可不再是bsites(1,b)的那个b(b的前n个与i=1~nn同步)，而是cut部分的<u>第几个点</u>

belongA就像是一个数组用来标记访问，对于cut的那一部分的每个格点标记为1

inversA的index是点的坐标，像是先获取到这个点的坐标，就能知道这个点是cut部分的第几个点.



别忘了还有一个 jpinsplit(nn) 用来标记A中属于B的格点.



```fortran
in subroutine adjustensemble(lambda):
do i=1,al:
if (spin(i)==spin1(i)) then
                 joinorsplit(i)=1 or 0
->
do i=1,al:
if (spin(cutsites(i))==spin1(cutsites(i))) then
                 joinorsplit(i)=1 or 0
```

这样看起来，关键点在于如何将切割写进 1 to al 这个范围中,  在corner这个代码里用的是al代指切割部分的点序号，比如切了的总共有4个点，那么al=4.

在 cornerless 的切割方式中，al = Lx*Ly/2 这是我们常见的. 如何理解？——和上面是一样的道理！

现在似乎理解了cutsites的作用了，它会将 1 to al 读取为cut的部分的第几个点，输出对应的点的坐标计数.

```fortran
in subroutine diagonalupdate():
->
do i=1,al
    copyspin(i)=spin(cutsites(i))
enddo

do i=1,al
     if (joinorsplit(i)==1) then
         spin(cutsites(i))=copyspin(i) 
     endif
enddo
```

果然，在cutsites的坐标转换下，我们在对角更新的时候，为了保留cut部分的spin状态，仍然可以只使用 1 to al 范围的数组复制

在Loopupdate缝裤腿的过程中，belongA这个数组可以用来替换掉原来用作区分的 1 to al, al+1 to nn.

```fortran
in subroutine totalloopupdate():

!第一条腿的裤腿缝合
do s1=al+1,nn 
    v1=frstspinop(s1)
    if (v1/=-1) then
        v2=lastspinop(s1)
        vertexlist(v2)=v1
        vertexlist(v1)=v2
    endif
 enddo
->
do s1=1,nn
    if(belongA(s1)==0) then
        v1=frstspinop(s1)
        if (v1/=-1) then
                v2=lastspinop(s1)
                vertexlist(v2)=v1
                vertexlist(v1)=v2
        endif
    endif
 enddo
 
 !第二条裤腿的初始化
 ->
 do i=1,nn
    if(belongA(i)==0) then
        frstspinop1(i)=-1  
        lastspinop1(i)=-1  
    endif
 enddo
 
 !第二个replica的自旋更新：
 
 !Abar中
if(i>al) then
    if (frstspinop1(i)/=-1) then
            if (vertexlist(frstspinop1(i))==-1) spin1(i)=-spin1(i)
    else
            if (ran()<0.5) spin1(i)=-spin1(i)
    endif
  ->
 if(belongA(i)==0) then     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (frstspinop1(i)/=-1) then
            if (vertexlist(frstspinop1(i))==-1) spin1(i)=-spin1(i)
    else
            if (ran()<0.5) spin1(i)=-spin1(i)
    endif
    
!A中但不属于B的部分
  if(i<=al .and. joinorsplit(i)==0) then 
    if (frstspinop1(i)/=-1) then
        if (vertexlist(frstspinop1(i))==-1) spin1(i)=-spin1(i)
    else
        if (ran()<0.5) spin1(i)=-spin1(i)
    endif
      
 ->
 if(joinorsplit(inversA(i))==0) then
    if (frstspinop1(i)/=-1) then
        if (vertexlist(frstspinop1(i))==-1) spin1(i)=-spin1(i)
    else
        if (ran()<0.5) spin1(i)=-spin1(i)
    endif
            
!A中所属B的部分的第二个replica应该与第一个replica的自旋保持一致:
->
do s1=1,al
     if(joinorsplit(s1)==1) then
          spin1(cutsites(s1))=spin(cutsites(s1))
     endif
 enddo
 
!B部分的自旋遍历算符序列依次更新
do i=0,mm-1
    op=opstring(i)
    if (mod(op,2)==1) then
       b=op/2
       s1=bsites(1,b)
       s2=bsites(2,b)
       if (belongA(s1)==1) then
           if(joinorsplit(inversA(s1))==1) then
               spin1(s1)=-spin1(s1)
           endif
       endif
       if (belongA(s2)==1) then   
           if(joinorsplit(inversA(s2))==1) then
               spin1(s2)=-spin1(s2)
           endif
       endif
    endif
enddo
```

reference: 
Dongxi Liu's code
