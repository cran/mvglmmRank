\name{mvglmmRank-package}
\alias{mvglmmRank-package}
\docType{package}
\title{
Multivariate Generalized Linear Mixed Models for Ranking Sports Teams
}
\description{
Maximum likelihood estimates are obtained via an EM algorithm with either a first-order or a fully exponential Laplace approximation.
}
\details{
\tabular{ll}{
Package: \tab mvglmmRank\cr
Type: \tab Package\cr
Version: \tab 1.2-2\cr
Date: \tab 2018-08-23\cr
License: \tab GPL-2\cr
}
See the help pages for mvglmmRank and game.pred
}
\author{
Andrew T. Karl and Jennifer Broatch

Maintainer: Andrew T. Karl <akarl@asu.edu>
}
\references{
Broatch, J.E. and Karl, A.T. (2017). Multivariate Generalized Linear Mixed Models for Joint Estimation of Sporting Outcomes. \emph{Italian Journal of Applied Statistics}. To appear. Also available from https://arxiv.org/abs/1710.05284.

Karl, A.T., Yang, Y. and Lohr, S. (2013). Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to 
    Educational Value-Added Assessments. \emph{Computational Statistics and Data Analysis},  59, 13-27.

Karl, A., Yang, Y. and Lohr, S. (2014) Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Karl, A.T. (2012). The Sensitivity of College Football Rankings to Several Modeling Choices, \emph{Journal of Quantitative Analysis in Sports}, 
    Volume 8, Issue 3, DOI 10.1515/1559-0410.1471
}
\keyword{ package }
\examples{
data(nfl2012)
mvglmmRank(nfl2012,method="PB0",first.order=TRUE,verbose=TRUE,max.iter.EM=1)
\donttest{
result <- mvglmmRank(nfl2012,method="PB0",first.order=TRUE,verbose=TRUE)
print(result)
game.pred(result,home="Denver Broncos",away="Green Bay Packers")
}
}
