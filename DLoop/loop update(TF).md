# Loop update

heat bath solution gives:  

$$
p_i = \frac{W_i}{\sum_i W_i}
$$

它表示从vertex的某个确定的入口进入后，流向下一个vertex $W_i$ 所选择的出口的概率，下一个vertex的权重比上所有从这个入口进入后从所有出口出能流向的vertex的权重之和.  

site的热浴解:   

$$
W_{off-diag} = gamma, W_{diag} = wxd
$$

$$
a_{off-diag} = ao = \frac{W_{off-diag} \times W_{off-diag}}{W_{off-diag}+wxd}
$$

$$
b_{off-diag} = bo = \frac{W_{off-diag} \times wxd}{W_{off-diag}+wxd}
$$

$$
a_{diag} = ad = \frac{wxd \times wxd}{W_{off-diag}+wxd}
$$

$$
b_{diag} = bd = \frac{W_{off-diag} \times wxd}{W_{off-diag}+wxd}
$$

```fortran
 wxo = gamma
 wxd = gamma !取任意权重，但最好不要太小，否则local更新插不进去
 ao = (wxo*wxo)/(wxo + wxd) !非对角site 反弹
 bo = (wxo*wxd)/(wxo + wxd) !非对角site 停止
 ad = (wxd*wxd)/(wxo + wxd) !对角site 反弹
 bd = (wxd*wxo)/(wxo + wxd) !对角site 停止
```


然后在对角更新处也有相应的修改:  
```fortran
 aprobg = wxd*beta*narray
 dprobg = 1.0/aprobg
```
因为对角更新里插入的是对角site算符，只需要用对角site的权重.

------

#### 在loop update中使用概率  
当遇到一个site算符，对其进行判断

```fortran
