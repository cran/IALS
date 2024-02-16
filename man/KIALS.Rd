\name{KIALS}
\alias{KIALS}
\title{
Estimating the Pair of Factor Numbers via Eigenvalue Ratios Corresponding to IALS
}
\description{
The function is to estimate the pair of factor numbers via eigenvalue ratios corresponding to IALS method. 
}
\usage{
KIALS(X, W1 = NULL, W2 = NULL, kmax, max_iter = 100, ep = 1e-06)
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
  \item{kmax}{
  The user-supplied maximum factor numbers. Here it means the upper bound of the number of row factors and column factors.
}
  \item{max_iter}{
  The maximum number of iterations for the algorithm, default is 100. See in \code{\link{IALS}}.
}
  \item{ep}{
  The stopping criterion in the iteration algorithm, default is \eqn{10^{-6} \times Tp_1 p_2}. See in \code{\link{IALS}}.
}
}
\details{
In detail, we first set \eqn{k_{\max}} is a predetermined upper bound for \eqn{k_1,k_2} and thus by \code{\link{IALS}} method, we can obtain the estimate of \eqn{\bold{F}_t}, denote as \eqn{\hat{\bold{F}}_t}, which is of dimension \eqn{k_{\max}\times k_{\max}}. Then the dimensions \eqn{k_1} and \eqn{k_2} are further determined as follows:
\deqn{\hat{k}_{1}=\arg\max_{j \leq k_{\max}}\frac{\lambda_{j}\left(\dfrac{1}{T}\sum_{t=1}^{T}\hat{\bold{F}}_t\hat{\bold{F}}_t^\top\right)}{\lambda_{j+1}\left(\frac{1}{T}\sum_{t=1}^{T}\hat{\bold{F}}_t\hat{\bold{F}}_t^\top\right)},}
\deqn{\hat{k}_{2}=\arg\max_{j \leq k_{\max}}\frac{\lambda_{j}\left(\dfrac{1}{T}\sum_{t=1}^{T}\hat{\bold{F}}_t^\top\hat{\bold{F}}_t\right)}{\lambda_{j+1}\left(\frac{1}{T}\sum_{t=1}^{T}\hat{\bold{F}}_t^\top\hat{\bold{F}}_t\right)}.}
}
\value{
\item{\eqn{k_1}}{The estimated row factor number.}
\item{\eqn{k_2}}{The estimated column factor number.}
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

kmax=8
K=KIALS(X, W1 = NULL, W2 = NULL, kmax, max_iter = 100, ep = 1e-06);K
}