# 1D Transverse field Ising model
FM Hamiltonian:  
$$H=-J\sum(\sigma^z_{i}\sigma^z_{j})-h\sum\sigma^x_{i}$$  
construct the element:  
$$H_{0,0}=I$$
$$H_{-1,g}=h(\sigma^+ + \sigma^-)$$
$$H_{0,g} = h$$
$$H_{1,b} = J(\sigma^z_i\sigma^z_j + 1)$$  
According detailed balance, we obtain:
$$P_d=\frac{M-n+1}{\beta(Nh+2JN_b)}$$
$$P_i=\frac{\beta(Nh+2JN_b)}{M-n}$$
insert a bond:
$$P_{bond}=\frac{2J N_b}{Nh+2JN_b}$$
energy:
$$E_{av}\sim\frac{\langle n\rangle}{-\beta}+J+h$$
