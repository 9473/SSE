非平衡方法测量 Rényi 纠缠熵

Definition of n-th Rényi entropy:


$$
S^n = \frac{1}{1-n}\ln{\frac{Z^n_A}{Z^n_\emptyset}}
$$


其中 $Z^n_{\emptyset}$ 表示的是 n个replica 的整个系统的配分函数，而 $Z^n$ 表示的 n个replica 的 subsystem A 的配分函数。  

非平衡的会引入一个参数 $\lambda$ . 即,  


$$
Z_A^n(\lambda) = \sum_{B \text{ in }A} \lambda^{N_B} (1-\lambda)^{N_A-N_B} Z_B^n
$$


我们设函数 $g(\lambda,N_B) = \lambda^{N_B} (1-\lambda)^{N_A-N_B}$ , 于是 $Z_A^n(\lambda) = \sum  g(\lambda,N_B) Z_B^n$.   B 是我们在 subsystem A 中挖的一个小部分. (这个说法可能有点不对，具体想象应该是 B是A与环境的glued geometry) 


$$
Z^n = e^{-\beta [-\frac{1}{\beta} \ln (g  + Z^n_B)]} = e^{-\beta W^n}
$$


  

  

由于引入的这么一个参数, $F(1)-F(0) = \int^1_0 f(x)dx$ , 熵的定义式子可以重新写成：  


$$
S^n =\frac{1}{1-n} \int^1_0 d\lambda \frac{\partial \ln Z^n_A(\lambda)}{\partial \lambda}
$$


这是因为  $Z^n_{\emptyset} = Z^n(\lambda = 0)$ 可以看成是 $Z^n_A(\lambda)$ 的一个特殊情况,  当 $\lambda = 1$ 时相当于 subsystem A 全部粘在 环境 $\bar{A}$ 上.  



仅仅是获得这样一个积分式子是不够的.   前文献中就对此讨论，因为QMC实际利用这个式子计算纠缠熵的操作有限，需要在 $\lambda = 0$ 和 $\lambda = 1$ 之间的细网格点上进行独立的平衡模拟，再对所有模拟结果的曲线进行数值积分。  

因此，更进一步地我们将 $\lambda$ 再进行细分.  $\lambda \to \lambda(t)$.  

  

这样熵的积分式子根据链式法则变成：  


$$
S^n =\frac{1}{1-n} \int^{t_f}_{t_i} dt \frac{d \lambda}{ dt} \frac{\partial \ln Z^n_A(\lambda)}{\partial \lambda}
$$

其中 $\lambda(t_f)=1, \lambda(t_i) = 0$.  

  

为了在有限的淬火时间里也能准确估计纠缠熵，我们还需要用到 Jarzynski 等式:  


$$
\Delta F = -\frac{1}{\beta} \ln \langle e^{-\beta \delta W}\rangle
$$


它的物理意义是，可以通过非平衡测量的系综平均值，通过以上关系，来提取在该路径中累积的自由能差 $\Delta F$.  由于 $ F = -\frac{1}{\beta} \ln Z \to \Delta F = -\frac{1}{\beta}\ln(\Delta Z)$ .  

因此, $\ln$ 内部之差变成 $\ln$ 之比, 利用熵的定义式子.  熵可以重新写成:  


$$
S^n = \frac{1}{1-n}\ln \langle e^{-\beta \delta W^n}\rangle
$$


蒙特卡洛实现的就是对括号内部的做系综平均.  

  



这其中还有一些trick.  

由于,  我们其实已经把 $\lambda$ 靠 $t$ 分成了很多片，我们假设这些片是均匀的，每一片宽度是 $\Delta$, 我们这里与文献保持一致，采用 $k$ .


$$
\frac{Z^n_A}{Z^n_\emptyset} = \frac{ Z^n(\lambda = 1)}{ Z^n(\lambda = 0)} = \frac{Z^n(\lambda = \Delta)}{Z^n(\lambda = 0)}\frac{Z^n(\lambda = 2\Delta)}{Z^n(\lambda = \Delta)} ...\frac{Z^n(\lambda = 1)}{Z^n(\lambda = 1-\Delta)}
$$

$$
= \prod_{k =1}^{K} \frac{Z^n(k\Delta)}{Z^n[(k-1)\Delta]}
$$


相当于 份数 累乘直到 $K = \frac{1}{\Delta}$.  因而,  


$$
S^n = \frac{1}{1-n}\ln \langle \prod_{k =1}^{K}e^{-\beta \delta W_k^n}\rangle
$$

文献中将 $\prod$ 放至 $\ln$ 前面变成 $\sum$, 但本质说的是一件事.  

  

我们之前已经做了 $Z^n = e^{-\beta [-\frac{1}{\beta} \ln (g  + Z^n_B)]} = e^{-\beta W^n}$ ,  


$$
W_k^n = -\frac{1}{\beta}\ln {[g(\lambda_k,N_B^k)  + Z^n_B]}
$$  
  
$\delta W_k^n = ?$. 


$$
\delta W_k^n = H(t+\delta t)-H(t) 
$$

$$
=-\frac{1}{\beta}\ln \frac{g(\lambda_k,N_B^k)}{g(\lambda_{k-1},N_B^{k-1})}
$$



 
$$
S^n = \frac{1}{1-n}\ln \langle \prod_{k =1}^{K} e^{-\beta [-\frac{1}{\beta}\ln \frac{g(\lambda_k,N_B^k)}{g(\lambda_{k-1},N_B^{k-1})}]}\rangle =
\frac{1}{1-n}\ln \prod_{k =1}^{K} \langle \frac{g(\lambda_k,N_B^k)}{g(\lambda_{k-1},N_B^{k-1})} \rangle
$$


因此，我们现在是对 $\frac{g(\lambda_k,N_B^k)}{g(\lambda_{k-1},N_B^{k-1})}$ 进行采样做平均.



(2024.3.15) 我现在仍然对 Jarzynski 等式和熵对 $t$ 的积分式子之间的联系有迷惑.
