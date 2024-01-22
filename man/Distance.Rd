\name{Distance}
\alias{Distance}
\title{
The distance between the spaces spanned by the column of two matrices.
}
\description{
Calculate the distance between spaces spanned by the column of two matrices. The distance is between 0 and 1. If the two spaces are the same, the distance is 0. if the two spaces are orthogonal, the distance is 1.
}
\usage{
Distance(Z1, Z2)
}
\arguments{
  \item{Z1}{
  Input a matrix with \eqn{p \times q_1}.
}
  \item{Z2}{
  Input another matrix with \eqn{p \times q_2}
}
}
\details{
Define \deqn{\mathcal{D}(\bold{Q}_1,\bold{Q}_2)=\left(1-\frac{1}{\max{(q_1,q_2)}}\text{Tr}\left(\bold{Q}_1\bold{Q}_1^\top\bold{Q}_2\bold{Q}_2^\top\right)\right)^{1/2}.}

By the definition of \eqn{\mathcal{D}(\bold{Q}_1,\bold{Q}_2)}, we can easily see that \eqn{0 \leq \mathcal{D}(\bold{Q}_1,\bold{Q}_2)\leq 1}, which measures the distance between the column spaces spanned by two orthogonal matrices \eqn{\bold{Q}_1} and \eqn{\bold{Q}_2}, i.e., \eqn{\text{span}(\bold{Q}_1)} and \eqn{\text{span}(\bold{Q}_2)}. In particular, \eqn{\text{span}(\bold{Q}_1)} and \eqn{\text{span}(\bold{Q}_2)} are the same when \eqn{\mathcal{D}(\bold{Q}_1,\bold{Q}_2)=0}, while  \eqn{\text{span}(\bold{Q}_1)} and \eqn{\text{span}(\bold{Q}_2)} are orthogonal when \eqn{\mathcal{D}(\bold{Q}_1,\bold{Q}_2)=1}. The Gram-Schmidt orthogonalization can be used to make \eqn{\bold{Q}_1} and \eqn{\bold{Q}_2} column-orthogonal matrices.
}
\value{
Output a number between 0 and 1.
}
\references{
He, Y., Zhao, R., & Zhou, W. X. (2023). Iterative Least Squares Algorithm for Large-dimensional Matrix Factor Model by Random Projection. <arXiv:2301.00360>.
}
\author{
Yong He, Ran Zhao, Wen-Xin Zhou.
}

\examples{
set.seed(1111)
A=matrix(rnorm(10),5,2)
B=matrix(rnorm(15),5,3)
Distance(A,B)
}
