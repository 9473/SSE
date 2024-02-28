# SSE
stochastic series expansion
link：https://tourmaline-sunscreen-c8f.notion.site/SSE-Stochastic-Series-Expans-QMC-209818ff11cc4bbbb8b0067c1830c130
more ditailed explanation

#### Pre work:

$$
H = J\sum_{\langle i,j\rangle} [S_i^xS_j^x+S_i^yS_j^y+S_i^zS_j^z]
$$

$$
{\color{Gray} H = -J\sum h, h = -[S_i^xS_j^x+S_i^yS_j^y+S_i^zS_j^z]} 
$$

$h$ makes no sense just write for compare  

$$
{\color{Purple} H = -J\sum_b^{N_b}H_b, H_b = H_{1,b}-H_{2,b} } 
$$
$$
H_{1,b} =-S_i^zS_j^z,\\-H_{2,b} =-[S_i^xS_j^x+S_i^yS_j^y]\to H_{2,b} = S_i^xS_j^x+S_i^yS_j^y
$$

$$
\to {\color{Purple} H_{2,b} = \frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j)}
$$

自旋交换作用  

$$
H_{1,b} = C-S_i^zS_j^z
$$

when $-S_i^zS_j^z |\uparrow \uparrow \rangle = -\frac{1}{4}|\uparrow \uparrow \rangle$ , for positive, $C = 1/4$  

$$
{\color{Purple}H_{1,b} = \frac{1}{4}-S_i^zS_j^z }
$$

伊辛相互作用  

$$
(\frac{1}{4}-S_i^zS_j^z)|\uparrow \uparrow \rangle = 0|\uparrow \uparrow \rangle,(\frac{1}{4}-S_i^zS_j^z)|\downarrow \uparrow \rangle = 1/2|\downarrow \uparrow \rangle,
$$
$$
(\frac{1}{4}-S_i^zS_j^z)|\uparrow \downarrow \rangle = 1/2|\uparrow \downarrow \rangle,(\frac{1}{4}-S_i^zS_j^z)|\downarrow \downarrow \rangle = 0|\downarrow \downarrow \rangle
$$

Only anti-paralleled term has contribution for Z  

$$
\frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j)|\uparrow \uparrow \rangle = 0,\frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j)|\downarrow \uparrow \rangle = 1/2|\uparrow \downarrow \rangle,
$$
$$
\frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j)|\uparrow \downarrow \rangle=1/2|\downarrow \uparrow \rangle ,\frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j)\downarrow \downarrow \rangle=0
$$  

Only anti-paralleled can act. So  

$$
{\color{Orange}H_{1,b}|\cdot*\rangle=1/2|\cdot*\rangle,H_{2,b}|\cdot*\rangle=1/2|*\cdot\rangle}
$$

------

#### Energy

$$
{\color{Purple}  Z = \sum_\alpha \sum _{n=0}^\infty \frac{(-\beta)^n}{n!}\bra{\alpha} H^n \ket{\alpha} }
$$

$$
\langle H \rangle  =\frac{1}{Z}\sum_\alpha \sum _{n=0}^\infty \frac{(-\beta)^n}{n!}\bra{\alpha} H^n H\ket{\alpha}
$$

$k = n-1 \to $  

$$
\langle H\rangle=\frac{1}{Z}\sum_\alpha \sum _{n=0}^\infty \frac{(-\beta)^k}{k!}\bra{\alpha} H^k H\ket{\alpha}
$$

$$
\to \langle H \rangle  =\frac{1}{Z}\sum_\alpha \sum _{n=0}^\infty \frac{(-\beta)^{(n-1)}}{(n-1)!}\bra{\alpha} H^n\ket{\alpha}
$$

$\frac{(-\beta)^{(n-1)}}{(n-1)!} = \frac{\frac{(-\beta)^{n}}{-\beta}}{(n-1)!\cdot n}\frac{n}{1} = \frac{n}{-\beta}\frac{(-\beta)^n}{n!} \to$  

$$
\langle H \rangle = \langle E \rangle =\frac{\sum_\alpha \sum_{n=0}^\infty \frac{n}{-\beta}\frac{(-\beta)^n}{n!} \bra{\alpha}H^n\ket{\alpha}}{ \sum_\alpha \sum _{n=0}^\infty \frac{(-\beta)^n}{n!}\bra{\alpha} H^n \ket{\alpha}}
$$  

$$
{\color{Orange} \langle E \rangle=\frac{\langle n \rangle}{-\beta} }
$$

$$
\langle H^2 \rangle = \frac{\langle n(n-1) \rangle}{\beta^2}
$$
$$
\langle H \rangle^2 =\frac{\langle n \rangle^2}{\beta^2}
$$
\langle H^2 \rangle-\langle H \rangle^2 = \frac{1}{\beta^2}(\langle n^2-n \rangle-\langle n \rangle^2 ) = \frac{1}{\beta^2}(\langle n^2 \rangle-\langle n \rangle^2 -\langle n \rangle)
$$

$$
C =\frac{\langle H^2 \rangle-\langle H \rangle^2}{T^2} =\beta^2 \frac{1}{\beta^2}(\langle n^2 \rangle-\langle n \rangle^2 -\langle n \rangle)
$$

$$
={\color{Orange}\langle n^2\rangle-\langle n\rangle^2-\langle n\rangle}
$$

average expansion order:  

$$
\langle n \rangle=-\beta\langle H\rangle
$$

$H = -J\sum_bH_b$, set $J =1$

$$
{\color{Purple} \langle n \rangle=\beta\langle \sum_b H_b \rangle = \beta N_b |E_b| },\ \ E_b =\langle H_b \rangle
$$

------

#### truncate at a maximum power M 截断

$$
{\color{Purple} \langle n \rangle=\beta\langle \sum_b H_b \rangle = \beta N_b |E_b| }
$$

$$
C ={\langle n^2 \rangle-\langle n \rangle^2 -\langle n \rangle}
$$

$E_b$可以看作每个键的平均能量，每个键的平均能量有限， $\langle n\rangle \propto \beta N$ ，并且，当低温下，比热 $C\to 0$,  $\langle n^2 \rangle-\langle n \rangle^2 -\langle n \rangle \to 0,\langle n^2 \rangle-\langle n \rangle^2\to \langle n \rangle$，你就会发现 $n$ 的分布展宽有限， $\sqrt{\langle n^2 \rangle-\langle n \rangle^2} \approx \sqrt{\langle n \rangle} \propto (\beta N)^{1/2}$ . So, we can imagine that, **$n$ has an expectation depends on size of system and far of position zero(position zero!!!) n has very small probability.**

因此，不仅是计算成本上的考虑，还考虑到$n$的实际分布，我们可以做截断。

------

#### Fill-in null operator 填充单位算符

Insert $(M-n)$  I operators.

本质上不同的构型切片会造成不同数量的非单位算符.

However, these operators do not occur in the evaluation of the partition function, and must be accounted for by dividing the final expression for Z by their contribution(考虑简并性)  

$$
\text{pick the placement  }:\binom{M}{n} = \frac{M!}{n!(M-n)!}
$$

$\frac{Z}{\frac{M!}{n!(M-n)!}}\to$  

$$
Z = \sum_{\alpha}\sum_{S_M} \frac{(-\beta)^n(M-n)!}{M!}\langle \alpha|\prod_{i=1}^M H_{b_i}|\alpha \rangle
$$

- where the sum over n is now implicitly included in the sampling of $S_M$: $\sum^M$

------

#### Weight

$$
Z = \sum_{i,n}W_i(n)
$$

我们发现，这样定义的了权重之后，配分函数相当于对各权重进行求和，那么配分函数的差异是完全可以体现在$W$上的  

截断到了$M$的话，  

$$
Z =\sum_{i,n}^MW_i(n), W_i(n) =\frac{(-\beta)^n(M-n)!}{M!}\langle \alpha|\prod_{i=1}^n H_{b_i}|\alpha \rangle
$$

对于内部算符数量、类型不同导致的配分函数的差异，由Weight可以体现出来  

同样的，你也可以定义  

$$
W_i(n) =\frac{(-\beta E_i)^n}{n!} \to W_i(n) =\frac{(-\beta)^n(M-n)!}{M!}
$$

差别不大

------

#### local update

$$
\pi(a)p(a \to b) A(a)=\pi(b)p(b \to a) A(b)
$$

对于[0,0] -> [1,b] 这个过程，总是会有 $1/N_b$ 的概率。在 Insert 这一情况中，这个过程作为正过程，因此 $\frac{[1,b] \to [0,0]}{[0,0] \to [1,b]} =N_b$ ；在 Delete 中这一过程作为逆过程，因此 $\frac{[0,0] \to [1,b]}{[1,b] \to [0,0]} = 1/N_b$ .

这一不均匀终究还是单位算符的操作和对角算符的操作不对成所导致的。  

$$
p = \frac{W_b}{W_a}=\frac{\frac{(-\beta)^{(n+1)}(M-n-1)!}{M!}(1/2)}{\frac{(-\beta)^n(M-n)!}{M!}} = \frac{\beta}{2}\frac{1}{(M-n)}
$$

$$
{\color{Blue}  \text{select}: N_b  \text{placement in}  \pi(a)}
$$
$$
{\color{Blue} \text{select prob}: 1/N_b }
$$

$$
{\color{Blue}  P = \frac{\pi(b)A(b)}{\pi(a)A(a)} =\frac{\beta}{2}\frac{1}{(M-n)} \cdot \frac{1}{\frac{1}{N_b}} =\frac{\beta}{2}\frac{N_b}{(M-n)}}
$$

###### Delete:

Non-I: n -> n-1  

$$
{\color{Gray} p = \frac{\frac{(-\beta)^{(n-1)}(M-n+1)!}{M!}(1/2)^{n-1}}{\frac{(-\beta)^n(M-n)!}{M!}(1/2)^n} = \frac{2(M-n+1)}{\beta} }
$$

$$
{\color{Red} P([1,b] \to [0,0]) = \frac{2(M-n+1)}{N_b \beta} }
$$
$$
=p'(\frac{p_{select}([0,0] \to [1,b])}{p_{select}([1,b] \to [0,0])}) =p'\cdot (1/N_b)
$$


