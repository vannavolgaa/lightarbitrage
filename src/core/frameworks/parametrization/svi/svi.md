# SVI / SSVI documentation 

Simple snipet of SVI and SSVI model. See paper Arbitrage-Free SVI volatiity surfaces (Gatheral & Jacquier, 2013)

## 1) Representation

The raw SVI model expresses the total variance smile as: 
$$w(k) = a + b \{p(k-m) + \sqrt{(k-m)^2 + \sigma^2}\}$$
The JW-SVI in terms of raw parameters is given by: 
- $v_t = \frac{a + b \{-pm + \sqrt{m^2 + \sigma^2}\}}{t}$, gives the ATM variance 
- $\psi_t = \frac{b}{2\sqrt{w_t}}(\frac{-m}{\sqrt{m^2 + \sigma^2}} + p)$, gives the ATM skew 
- $p_t = \frac{b}{\sqrt{w_t}}(1 - p)$, the slope of the left wing (put)
- $c_t = \frac{b}{\sqrt{w_t}}(1 + p)$, the slope of the right wing (call)
- $ \tilde{v}_t = \frac{a + b\sigma \sqrt{1-p^2}}{t}$, the minimum implied variance 

Where $w_t := tv_t$. Inversely, the raw parameters can be obtained from the JW-SVI representation: 
- $b = \frac{\sqrt{w_t}}{2}(c_t + p_t)$
- $p = 1 - \frac{p_t\sqrt{w_t}}{b}$
- $a = t\tilde{v}_t - b \sigma \sqrt{(1-p^2)}$
- $m = \frac{(v_t - \tilde{v}_t)t}{b\{-p+sign(\alpha)\sqrt{(1+\alpha^2)}-\alpha*(1-p^2)\}}$
- $\sigma = \alpha m$ or $\sigma = \frac{w_t - a}{b}$ if $m=0$


The Surface SVI expresses the total variance surface as: 
$$w(k, \theta_t) = \frac{\theta_t}{2}\{1 + p\phi(\theta_t)k + \sqrt{(\phi(\theta_t)k +p)^2 + (1-p^2)}\}$$
where $\phi(\theta_t) = \nu \theta_t^\gamma$ is the power law parametrization. $\theta_t = w(0)_t$ is the ATM total variance at maturity $t$. 
Under this framework, the ATM volatility skew is given by: 
$$\partial_{k} \sigma_{BS}(k,t) \mid_{k=0} = \frac{p\sqrt{\theta_t}}{2\sqrt{t}}\phi(\theta_t)$$
The JW-SVI in term of SSVI is given by: 
- $v_t = \theta_t/t$
- $\psi_t = \frac{1}{2}p\sqrt{\theta_t}\phi(\theta_t)$
- $p_t = \frac{1}{2}\sqrt{\theta_t}\phi(\theta_t)(1-p)$
- $c_t = \frac{1}{2}\sqrt{\theta_t}\phi(\theta_t)(1+p)$
- $ \tilde{v}_t = v_t*(1-p^2)$

## 2) Conditions for an arbitrage free surface 
### a) Arbitrage free volatility surface 
A volatility surface is free of static arbitrage
if and only if the following conditions are satisfied:
- it is free of calendar spread arbitrage: A volatility surface $w$ is free of calendar spread if $\partial_t w(k,t) \geq 0$
- each time slice is free of butterfly arbitrage: A slice is free of butterfly arbitrage if and only if $g(k) \geq 0$ for all $k \in \mathbb{R}$ and if $\lim_{k \rightarrow 0} d_+(k) = -\inf$ (see Breeden and Litzenberger, 1978 for the first condition). The second condition comes from the fact that a call price must converges to 0 as k converges to infinity. 

Where $g(k) = (1 - \frac{kw'(k)}{2w(k)})^2 - \frac{w'(k)^2}{4}(\frac{1}{w(k)} + \frac{1}{4}) + \frac{w''(k)}{2}$


### b) Static arbitrage free Surface SVI 
The SSVI volatility surface is free of butterfly arbitrage if the following conditions are satisfied for all $\theta \geq 0$: 
- $\theta \phi(\theta) (1 + \mid p \mid) < 4$
- $\theta \phi(\theta)^2 (1 + \mid p \mid) \leq 4$ 

The SSVI volatility surface is free of calendar spread arbitrage if and only if: 
- $\partial_t \theta_t \geq 0$ for all $t \geq 0$
- $0 \leq \partial_{\theta}\phi(\theta) \leq \frac{1}{p^2}(1+\sqrt{1-p^2})\phi(\theta)$ for all $\theta \geq 0$

### c) Butterfly arbitrage free SVI smile
An JW-SVI smile is free of butterfly arbitrage if the two following conditions are satisfied for $t>0$: 
- $\sqrt{w_t}\max(p_t, c_t) < 2$
- $(p_t + c_t)\max(p_t, c_t) \leq 2$

### d) Calendar spread arbitrage free between two SVI slices 
Recall that raw SVI is: 
$$w(k) = a + b \{p*(k-m) + \sqrt{(k-m)^2 + \sigma^2}\}$$
setting $u = k-m$ this leads to 
$$w'(k) = b \{p + \frac{u}{\sqrt{u^2 + \sigma^2}}\}$$
First note that $\lim_{\mid u \mid \rightarrow \infty} \frac{1}{\sqrt{1 + \frac{\sigma^2}{u^2}}} = 1$. Thus, for the right slope ($u>0$) we have: 
$$w_+'(k) = b \{p + \frac{1}{\sqrt{1 +  \frac{\sigma^2}{u^2}}}\}$$
which leads to the asymptotic right slope: 
$$\tilde{w}_+'(k) = \lim_{k \rightarrow \infty} w_+'(k) = b(p+1)$$
Similary, for the left slope ($u<0$) we get: 
$$w_-'(k) = b \{p - \frac{1}{\sqrt{1 +  \frac{\sigma^2}{u^2}}}\}$$
$$\tilde{w}_-'(k) = \lim_{k \rightarrow -\infty} w_-'(k) = b(p-1)$$
One can obtain the asymtotic bounds $\tilde{k}_+$ and $\tilde{k}_-$ from which the raw SVI behaves close to the asymtotically. The problem is: 
$$\mid w_{sign(u)}'(k) -  \tilde{w}_{sign(u)}'(k)\mid \leq \epsilon \tilde{w}_{sign(u)}'(k)$$
We get: 
$$\tilde{k}_{sign(u)} = m + sign(u) \sigma \sqrt{\frac{\{1 - \epsilon(1 + sign(u) p)\}^2}{1 - \{1 - \epsilon(1 + sign(u) p)\}^2}}$$
For two SVI slices $w_1(k)$ and $w_2(k)$, with $t1 < t2$, to not display calendar spread arbitrage the following conditions must hold: 
- $t_1 v_1 \leq t_2 v_2$ : condition on ATM total variance (similar to SSVI representation)
- $t_1 \tilde{v}_1 \leq t_2 \tilde{v}_2$ : condition on minimum total variance 
- $t_1 v_1 p_1 \leq t_2 v_2 p_2$ : asymptotic condition on left slope 
- $t_1 v_1 c_1 \leq t_2 v_2 c_2$ : asymptotic condition on right slope 

If those conditions holds, one last numerical check must be done on the grid $ G = [\tilde{k}_-, \tilde{k}_+]$. The following conditions must hold for all $k_i \in G$: 
$$w_1(k_i) \leq w_2(k_i) $$


## 3) Expressing the local volatility with SVI 
Gatheral shows that the local variance can be expressed in terms of total variance by modyfing the initial black scholes definition in terms of total variance $w(k)$ and log strike $k = \log{\frac{K}{F_t}}$ such that the local variance $\sigma(k)$ is given by: 
$$\sigma(k) = \frac{\partial_t w(k)}{g(k)}$$
Where $g(k)$ is the risk neutral density given in 2.a. $g(k)$ can be expressed in term of raw SVI parameters, recall that:
$$g(k) = (1 - \frac{kw'(k)}{2w(k)})^2 - \frac{w'(k)^2}{4}(\frac{1}{w(k)} + \frac{1}{4}) + \frac{w''(k)}{2}$$
where: 

$w(k) = a + b \{p(k-m) + \sqrt{(k-m)^2 + \sigma^2}\}$

$w'(k) = b\{p + \frac{k-m}{\sqrt{(k-m)^2 + \sigma^2}}\}$

$w''(k) = \frac{b\sigma^2}{\{(k-m)^2 + \sigma^2\}^{3/2}}$

The derivative of the total variance with respect to $t$ comes down to compute the derivatives of the raw SVI parameters expressed in terms of JW-SVI parameters, which are time $t$ dependant. Recall the total variance in terms of raw SVI: 
$$w(k) = a + b \{p(k-m) + \sqrt{(k-m)^2 + \sigma^2}\}$$
We can express the first derivative of $w(k)$ with respect to time such that: 
$$\partial_t w(k) = (\partial_t a) + p \{u (\partial_t b) + b (\partial_t u) \} + (\partial_t b) \sqrt{u^2 + \sigma^2} + \frac{b}{2\sqrt{u^2+\sigma^2}} \{2u(\partial_t u) + 2\sigma (\partial_t \sigma)\}$$
where: 

$u = k - m$ and $\partial_t u = -\partial_t m$

$\partial_t b = \frac{\sqrt{v_t}}{4\sqrt{t}}(c_t + p_t)$

$\partial_t p = 0$, since $p = 1 - \frac{p_t\sqrt{w_t}}{b} = 1 - \frac{2p_t}{p_t + c_t}$ 

$\partial_t m = \frac{m(b - t (\partial_t b))}{tb}$

if $m = 0$: $\partial_t a = \frac{a}{t}$ and $\partial_t \sigma = \frac{b(v_t - (\partial_t a))+(\partial_t b)(w_t-a)}{b^2}$ 

else: $\partial_t \sigma = \alpha (\partial_t m)$ and $\partial_t a = \tilde{v}_t - \sqrt{1-p^2}\{\sigma(\partial_t b) + b(\partial_t \sigma)\}$

## 4) Reduced SVI
Reduced SVI is a modified slice version of the SVI express from three parameters $v_t$, $\eta$ and $\rho$. This comes from the reprsentation of the SSVI into JW-SVI parameters. One can note that by fixing $\eta = \sqrt{\theta_t}\phi(\theta_t) >= 0$ we get the following version of the JW-SVI representation: 
- $\psi_t = \frac{1}{2}\eta\rho$
- $p_t = \frac{1}{2}\eta(1-\rho)$
- $c_t = \frac{1}{2}\eta(1+\rho)$
- $ \tilde{v}_t = v_t*(1-\rho^2)$

This allows to reduce the slice parametrization to 3 parameters from the initial 5, which could help for calibration and stability of the model without loosing fexibility. 