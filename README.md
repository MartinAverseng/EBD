# EBD
Fast convolution with radial kernels in 2D, EBD toolbox (matlab)




This toolbox implements the algorithm described in "Discrete convolution in $\mathbb{R}^2$ 
with radial kernels using non-uniform fast Fourier transform with non-equispaced frequencies", 
written by Martin Averseng, and submitted to the journal Numerical Algorithms in 2018. 

To test it, you can directly run Demo.m 
Here follows a more detailed description of the algorithm and a tutorial.

The method is designed to compute fast approximations of vectors $q$ which entries are given by 
$$q_k = \sum_{l=1}^{N_y} G(X_k - Y_l) f_l, \quad k = 1, \cdots, N_x$$
where $G$ is a radial function i.e. $G(x) = g(|x|)$ for some function $g$, $X$ and $Y$ 
are two clouds of $N_x$ and $N_y$ points in $\mathbb{R}^2$ and $f$ is a complex 
vector. 

1°) Description of the algorithm

The method first decomposes $G$ in finite Bessel series 
$$G(r) = \sum_{p = 1}^P \alpha_p J_0(\rho_p r) $$
where $J_0$ is the Bessel function of first kind, $\rho$ is the sequence of its
positive zeros, and $\alpha$ are called the EBD (Efficient Bessel Decomposition) of $G$. 
The coefficients are chosen as the minimizers of the $H^1_0$ error in this approximation on 
a ring $r_{min} < r < r_{max}$ where rmin is a cutoff parameter and rmax is the greatest distance 
between two points $X_k$ and $Y_l$. The method takes the paramter $a:= \frac{r_{min}}{r_{max}}
as an input. 

Then, each $J_0(\rho_p r)$ is approximated by 
$$J_0(\rho_p |x|) = \frac{1}{M_p}\sum_{m = 1}^{M_p} e^{i \rho_p \xi_{m}^p \cdot x}$$
where $\xi_m^p = e^{i\frac{2m\pi}{M_p}}$. This is the trapezoidal rule applied to the 
identity 
$$J_0(|x|) = \int_{\partial B} e^{i x \cdot \xi} d\xi$$
where the integration takes place on the boundary of the unit ball $B$ in $\mathbb{R}^2$. 

Combining these two steps, we obtain an approximation for $G$ of the form 
$$G(x) \approx \sum_{\nu = 1}^{N_\xi} \hat{\omega}_\nu e^{i \xi_\nu \cdot x}$$
valid for $|x| > r_{min}$. If this is replaced in the expression of $q_k$, we see that 
$q$ can be approximated by non-uniform Fourier transform for any vector $f$. 
The interactions $|X_k - Y_l| < rmin$ where the approximation is not valid
are corrected by a sparse matrix product. 


The code is in Matlab language. The implementation of the NUFFT is borrowed from
Leslie Greengard, June-Yub Lee and Zydrunas Gimbutas (see license file in the libGgNufft2D folder).  
The ideas come from a similar method in 3D called Sparse Cardinal Sine Decomposition, 
developped by François Alouges and Matthieu Aussal, also published in Numerical Algorithms. 

2°) Tutorial

- Create the arrays X and Y of size Nxx2 and Nyx2 (points in R^2). 

- Create a kernel by calling 
> G = Kernel(fun,der)
where fun is an anonymous function of your choice and der is an anonymous function 
repesenting the derivative of fun. For example, 
> G = Kernel(@(x)(1./x),@(x)(-1./x.^2))
Note the functions must accept arrays as input.
You can also use particular kernels like
> G = LogKernel  % represents G(x) = log(x)
> G = ThinPlate(a,b) % represents G(x) = a*x^2*log(b*x)
(and others, see folder Kernels) for which the method has been optimized.

- Define the tolerance in the error of approximation. The method guarantees that 
the Bessel decomposition of $G$ in the ring is accurate at the tolerance level. 
This implies that the max error on an entry of $q_k$ is tol*norm(q,1).  

- Define the parameter $a$ (the ratio between $r_{min}$ and $r_{max}$.)
$a$ is roughly the proportion of interactions that will be computed exactly. 
There is an optimal value of $a$ for which the evaluation of the convolution 
is the fastest, but it is not possible to know it in advance. However, when $X$ and
$Y$ are uniformly distributed on a disk, the optimal $a$ is of the order
$\frac{1}{\sqrt{\sqrt{N_x N_y}}}$, and if they are uniformly distributed on
a curve, the optimal $a$ is of the order $\frac{1}{\sqrt{N_x N_y}^{2/3}}$. 
If $a$ is large, a lot of interactions are computed exactly while the Bessel decomposition 
will have only a few terms. If $a$ is small, the opposite will happen. 

- You can now call
>[onlineEBD, rq, loc] = offlineEBD(G,X,Y,a,tol)
onlineEBD is a callable object. If f is a vector with length equal to size(Y,1)
> onlineEBD(f) 
returns the approximation of the convolution. rq is a RadialQuadrature object, 
you can see how the approximation goes by calling
> rq.show()
and also check all its propertis (including coefficients, frenquencies used, accuracy, ...). 
loc is the local correction matrix (used inside onlineEBD). It can be used for preconditioning purposes. 




