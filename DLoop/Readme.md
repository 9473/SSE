# Directed Loop Algorithm

矩阵元的权重为  

$$
diag,bond,for \ anti: \Delta /2 + \varepsilon + h_b
$$

$$
off-diag,bond,anti:1/2
$$

$$
diag, bond,allup:\varepsilon + 2h_b
$$

$$
diag,bond,alldn:\varepsilon
$$

$$
site:gamma
$$

无论是diag site还是off-diag site产生的权重都是gamma，这种处理见横场Ising. 其中 $h_b=h_z/2,gamma=h_x/2$ in 1D.  

总哈密顿量形式为：  

$$
H = \sum (\Delta S^zS^z + \frac{1}{2}(S^+S^- + S^-S^+) - h_z \sum S^z - h_x\sum S^x 
$$


## local update(diagonal update)
涉及多种算符类型的对角更新，在这里有bond算符和site算符（site算符的引入见横场Ising）  


```fortran
narray = nb + nn
```

这一段的目的是将键算符和格点算符总共能插的位置排列（顺序上谁在前谁在后无所谓）

掷一个随机数

```fortran
s=min(int(ran()*narray)+1,narray)
```

如果随机数掷到bond算符位置范围，说明选定了bond算符类型，按  

$$
\frac{N_{array}\beta \langle H_{diag}\rangle }{M-n}
$$  

的概率插入相应的算符

```fortran
s=min(int(ran()*narray)+1,narray)
          if (s <= nb) then
             b = s !此时s就是要插入键的位置，设为b
             
             if (spin(bsites(1,b))/=spin(bsites(2,b))) then !在自旋相反的地方插入对角bond
                if (aprob>=dfloat(mm-nh).or.aprob>=ran()*(mm-nh)) then
                !here aprob=w2*beta*narray, w2=z/2+h+e
                   opstring(i)=4*b+1
                   nh=nh+1
                endif
             else
                  !如果spin一样，那要判断朝上朝下
                  if(spin(bsites(1,b))==-1) then
                     if (aprob1>=dfloat(mm-nh).or.aprob1>=ran()*(mm-nh)) then
                     !here aprob1=w3*beta*narray, w3=e
                        opstring(i)=4*b+1
                        nh=nh+1
                     endif
                  else
                     if (aprob2>=dfloat(mm-nh).or.aprob2>=ran()*(mm-nh)) then
                     !here aprob2=w4*beta*narray, w4=e+2*h
                        opstring(i)=4*b+1
                        nh=nh+1
                     endif
                  endif
               endif
               
           !all this diag bond: op=4*b+1
```

或者，掷到site算符：  

```fortran
      else 	!! if nb<s
           if (aprobg>=ran()*dfloat(mm-nh)) then
           !here aprobg = wx*beta*narray, wx = gamma
           g = s - nb !此时s-nb才是要插入的site位置
           opstring(i) = 4*g
           nh = nh+1
           endif
        endif 
       !diag site: op=4*g+0
```

如果检测到diag bond算符，按1/aprob的概率删去

```fortran
elseif(mod(op,4)==1) then !op=4*b+1
            b=op/4
            if (spin(bsites(1,b))/=spin(bsites(2,b))) then
               p=dprob*(mm-nh+1)
                  if ((p>=1.d0).or.(p>=ran()))then
                  opstring(i)=0
                  nh=nh-1
                  endif
            else!如果spin一样，那要判断朝上朝下
               if(spin(bsites(1,b))==-1) then
               p=dprob1*(mm-nh+1)
                  if ((p>=1.d0).or.(p>=ran())) then
                  opstring(i)=0
                  nh=nh-1
                  endif
               else
               p=dprob2*(mm-nh+1)
                  if ((p>=1.d0).or.(p>=ran())) then
                  opstring(i)=0
                  nh=nh-1
                  endif
               endif
            endif
```

又或者，检测到diag site算符，

```fortran
elseif((mod(op,4)==0)) then  !op=4*g
            p = dprobg*(mm-nh+1)
            if(p>=1.d0 .or. p>=ran()) then
                opstring(i) = 0
                nh=nh-1
            endif
```


最后，检测到非对角算符，相应的自旋进行翻转

```fortran
elseif((mod(op,4)==3)) then !op=4*g+3
!this is off diag site
    g=op/4
    spin(bsites(1,g))=-spin(bsites(1,g)) !only change spin(i)

else !mod=2
!this is off-diag bond
    b=op/4
    spin(bsites(1,b))=-spin(bsites(1,b))
    spin(bsites(2,b))=-spin(bsites(2,b))
endif
```

在这里，opstring的计数方式为 $H(a,b)$ :  

$$
H_{1,b} = diag,bond
$$

$$
H_{2,b} = off-diag,bond
$$

$$
H_{0,g} = diag,site
$$

$$
H_{3,g} = off-diag,site
$$  

这样的话, 可以直接通过 mod(op,4) 来进行判断算符类型:  

$$
opstring(i) = op, op = 4*b+a, a = 0,1,2,3
$$

$$
a=1,op=4*b+1,b=op/4
$$

$$
a=2,op=4*b+2,b=op/4
$$

$$
a=3,op=4*g+3,g=op/4
$$  

------

