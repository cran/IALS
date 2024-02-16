\name{IALS}
\alias{IALS}
\title{
Iterative Alternating Least Square Estimation for Large-dimensional Matrix Factor Model
}
\description{
This function is designed to fit the matrix factor model using the Iterative Least Squares (IALS) method, rather than Principal Component Analysis (PCA)-based methods. In detail, in the first step, we propose to estimate the latent factor matrices by projecting the matrix observations with two deterministic weight matrices, chosen to diversify away the idiosyncratic components. In the second step, we update the row/column loading matrices by minimizing the squared loss function under the identifiability condition. The estimators of the loading matrices are then treated as the new weight matrices, and the algorithm iteratively performs these two steps until a convergence criterion is reached.
}
\usage{
IALS(X, W1 = NULL, W2 = NULL, m1, m2, max_iter = 100, ep = 1e-06)
}
\arguments{
  \item{X}{
  Input an array with \eqn{T \times p_1 \times p_2}, where \eqn{T} is the sample size, \eqn{p_1} is the the row dimension of each matrix observation and \eqn{p_2} is the the column dimension of each matrix observation.
}
  \item{W1}{
  The initial value for the row factor loading matrix. The default is NULL, with an initial estimate chosen from \eqn{\alpha}-PCA if not provided.
}
  \item{W2}{
  The initial value for the column factor loading matrix. The default is NULL, with an initial estimate chosen from \eqn{\alpha}-PCA if not provided.
}
  \item{m1}{
  A positive integer indicating the row factor number.
}
  \item{m2}{
  A positive integer indicating the column factor number.
}
  \item{max_iter}{
  The maximum number of iterations for the algorithm, default is 100.
}
  \item{ep}{
  The stopping criterion in the iteration algorithm, default is \eqn{10^{-6} \times Tp_1 p_2}.
}
}
\details{
Assume we have two weight matrices \eqn{\bold{W}_i} of dimension \eqn{p_{i} \times m_{i}} for \eqn{i=1,2}, as substitutes for \eqn{\bold{R}} and \eqn{\bold{C}} respectively. Then it is straightforward to estimate \eqn{\bold{F}_t} simply by
\deqn{\hat{\bold{F}}_{t}=\frac{1}{p_{1}p_{2}}\bold{W}_1^\top\bold{X}_t\bold{W}_2.}

Given \eqn{\hat{\bold{F}}_t} and \eqn{\bold{W}_1}, we can derive that 
\deqn{\hat{\bold{R}}=\sqrt{p_{1}}\left(\sum_{t=1}^{T}\bold{X}_t\bold{W}_2\hat{\bold{F}}_t^\top\right)\left[\left(\sum_{t=1}^{T}\hat{\bold{F}}_t\bold{W}_2^\top\bold{X}_t^\top\right)\left(\sum_{t=1}^{T}\bold{X}_t\bold{W}_2\hat{\bold{F}}_t^\top\right)\right]^{-1/2}.}

Similarly, we get the following estimator of the column factor loading matrix
\deqn{\hat{\bold{C}}=\sqrt{p_{2}}\left(\sum_{t=1}^{T}\bold{X}_t^\top\hat{\bold{R}}\hat{\bold{F}}_t\right)\left[\left(\sum_{t=1}^{T}\hat{\bold{F}}_t^\top\hat{\bold{R}}^\top\bold{X}_t\right)\left(\sum_{t=1}^{T}\bold{X}_t^\top\hat{\bold{R}}\hat{\bold{F}}_t\right)\right]^{-1/2}.}

Afterwards, we sequentially update \eqn{\bold{F}}, \eqn{\bold{R}} and \eqn{\bold{C}}. In simulation, the iterative procedure is terminated either when a pre-specified maximum iteration number (\eqn{\text{maxiter}=100}) is reached or  when \deqn{\sum_{t=1}^{T}\|\hat{\bold{S}}^{(s+1)}_t-\hat{\bold{S}}^{(s)}_t\|_{F} \leq \epsilon \cdot Tp_1p_2,}
where \eqn{\hat{\bold{S}}^{(s)}_t} is the common component estimated at the \eqn{s}-th step, \eqn{\epsilon} is a small constant (\eqn{10^{-6}}) given in advance.
}
\value{
The return value is a list. In this list, it contains the following:
\item{R}{The estimated row loading matrix of dimension \eqn{p_1\times m_1}, satisfying \eqn{\bold{R}^\top\bold{R}=p_1\bold{I}_{m_1}}.}
\item{C}{The estimated column loading matrix of dimension \eqn{p_2\times m_2}, satisfying \eqn{\bold{C}^\top\bold{C}=p_2\bold{I}_{m_2}}.}
\item{F}{The estimated factor matrix of dimension \eqn{T \times m_1\times m_2}.}
\item{iter}{The number of iterations when the stopping criterion is met.}
}
\references{
He, Y., Zhao, R., & Zhou, W. X. (2023). Iterative Alternating Least Square Estimation for Large-dimensional Matrix Factor Model. <arXiv:2301.00360>.
}
\author{
Yong He, Ran Zhao, Wen-Xin Zhou.
}

\examples{
set.seed(11111)
T=20;p1=20;p2=20
k1=3;k2=3

R=matrix(runif(p1*k1,min=-1,max=1),p1,k1)
C=matrix(runif(p2*k2,min=-1,max=1),p2,k2)

X=E=array(0,c(T,p1,p2))
F=array(0,c(T,k1,k2))

for(t in 1:T){
  F[t,,]=matrix(rnorm(k1*k2),k1,k2)
  E[t,,]=matrix(rnorm(p1*p2),p1,p2)
}
for(t in 1:T){
X[t,,]=R\%*\%F[t,,]\%*\%t(C)+E[t,,]
}

#Estimating the matrix factor model using the default initial values
fit1 = IALS(X, W1 = NULL, W2 = NULL,k1, k2, max_iter = 100, ep = 1e-06) 
Distance(fit1$R,R);Distance(fit1$C,C)

#Estimating the matrix factor model using one-step iteration
fit2 = IALS(X, W1 = NULL , W2 = NULL, k1, k2, max_iter = 1, ep = 1e-06) 
Distance(fit2$R,R);Distance(fit2$C,C)
}
