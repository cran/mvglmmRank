\name{mvglmmRank}
\alias{mvglmmRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mvglmmRank
}
\description{
This function fits a (multivariate) generalized linear mixed model to team scores and/or win/loss indicators.
}
\usage{
mvglmmRank(game.data, method = "PB0", first.order = FALSE, 
home.field = TRUE, max.iter.EM = 1000, tol1 = 1e-04, 
tol2 = 1e-04, tolFE = 0, tol.n = 1e-07, verbose = TRUE, OT.flag = FALSE, Hessian = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{game.data}{
a data frame that contains a column "home" of team names, a column "away" of team names,
a column "home.response" containing the scores (or other response) of the "home" teams, a column "away.response"
containing the scores (or other response) of the "away" teams, (optionally) a column "binary.response" that contains a column of binary responses (0's and 1's), and (optionally) a column "neutral.site" which 
takes the value 1 for neutral site games and 0 otherwise. 
NOTE: If game.data does not contain a "binary.response" column, then an indicator will be created
for whether the home team won.
NOTE: For neutral site games, randomly assign the teams as "home" or "away". As noted below, the data frame may optionally contain a column,
OT, which indicates how many overtime periods were played. 
NOTE: the game.data$OT column should not contain missing data. If there was no overtime, specify "none" or 0.
}
  \item{method}{
a character (remember to use quotation marks!). Choices are "N", "P0", "P1", "B",
"NB", "PB0", "PB1", or "NB.mov". "N" indicates the scores are fit with a normal distribution with intra-game
correlation between the home and away teams accounted for in an unstructured 2x2 error covariance matrix.
"P" indicates the scores are fit with a Poisson distribution. "B" indicates the home win/loss indicators are
fit using a binary distribution with a probit link. The presence of a "1" with a "P" indicates
potential intra-game correlation is modeled with an additional game-level random effect. A "0"
indicates no such random effects are included. "NB.mov" fits the margin of victory of the "home" team (under an assumed normal distribution) 
jointly with the binary win/loss indicators.
}
  \item{first.order}{
logical. TRUE requests that only a first order Laplace approximation be used, FALSE
requests a fully exponential Laplace approximation. See the references.
}
  \item{home.field}{
logical. TRUE requests that seperate home and away mean scores be modeled (along with a mean 
neutral site score, if applicable) along with a single home field effect in the binary model. 
FALSE requests only a single mean be calculated for the scores, and no fixed effects are fit for the
binary win/loss indicators.
}
  \item{max.iter.EM}{
a number giving the maximum number of EM iterations.
}
  \item{tol1}{
refers to the maximum relative change in parameters between iterations. This is the convergence criterion for the 
first order Laplace approximation. The first order Laplace approximation runs 
until tol1 signals, at which point the fully exponential corrections for the random effects vector begin
}
  \item{tol2}{
The fully exponential iterations run until the maximum relative change in model paramters is less than tol2. N/A when first.order==TRUE. 
}
  \item{tolFE}{
intermediate convergence criterion for fully exponential approximations. The algorithm runs with the fully exponential corrections only to the random
effects vector until tolFE signals (maximum relative change in parameters). 
After this, the fully exponential corrections for both the random efffects vector and the random effects covariance matrix are calculated 
}
  \item{tol.n}{
convergence tolerance for EM algorithm with method="N". Convergence is declared when \eqn{(l_k-l_{k-1})/l_k < tol.n}, where \eqn{l_k} is the log-likelihood 
at iteration \eqn{k}.
}
  \item{verbose}{
logical. If TRUE, model information will be printed after each iteration.
}
  \item{OT.flag}{
logical. If TRUE, then there should be a continuous column OT in game.data that indicates how many overtime periods there were for each game. The information will not be used for the binary models. NOTE: the game.data$OT column should not contain  missing data. If there was no overtime, specify 0.
}
\item{Hessian}{
logical. If TRUE, the Hessian of the model parameters is calculated via a central difference approximation.
}
}
\details{
Setting first.order=TRUE will yield the first order Laplace approximation. A partial fully exponential Laplace approximation can be obtained by setting tol1 > tol2 and tolFE=0. This
will apply fully exponential corrections to the vector of team ratings (the EBLUPs), but not to the covariance matrix of this vector. Karl, Yang, and Lohr (2014) show that this approach
produces a large portion of the benefit of the fully exponential Laplace approximation in only a fraction of the time. Using the default tolerances of mvglmmRank leads to this behavior.

To summarize, the models (except for method="N") run with the first order Laplace approximation until the relative change between parameteres is <= tol1. If first.order=TRUE, the program stops.
Otherwise, the program continues with the Laplace approximation, applying fully exponential corrections to the random effects vector until the maximum of the relative parameter
changes is <= tolFE. At this point, the program continues using the complete fully exponential Laplace approximation (corrections to both the random effects vector and its covariance matrix) 
until the maximum relative parameter change is <= tol2. If tolFE < tol2, then the program will finish without applying fully exponential corrections to the random effects covariance matrix.

method="PB1" is the least scalable, as the memory and computational requirements for this model are at least O((teams+number of games)^2). In the example data included with the package, the
NCAA basketball data is slow with the fully exponential approximation and method="PB1".
}
\value{
mvglmmRank returns an object of class \code{mvglmmRank}



An object of class \code{mvglmmRank} is a list containing the following components:

\item{n.ratings.offense}{The vector of offensive ratings from the normal model, or \code{NULL} if the normal model was not fit.}
\item{n.ratings.defense}{The vector of defensive ratings from the normal model, or \code{NULL} if the normal model was not fit.}
\item{p.ratings.offense}{The vector of offensive ratings from the Poisson model, or \code{NULL} if the Poisson model was not fit.}
\item{p.ratings.defense}{The vector of defensive ratings from the Poisson model, or \code{NULL} if the Poisson model was not fit.}
\item{b.offense}{The vector of win-propensity ratings from the binary model, or \code{NULL} if the binary model was not fit.}
\item{n.mean}{Mean scores from the normal model.}
\item{p.mean}{Mean scores from the Poisson model.}
\item{b.mean}{Home field effect from the binary model.}
\item{G}{Single block of random effects covariance matrix.}
\item{G.cor}{Correlation matrix corresponding to covariance matrix \code{G}.}
\item{R}{Error covariance matrix for normal model, or NULL if normal model not used.}
\item{R.cor}{Error correlation matrix for normal model, or NULL if normal model not used.}
\item{home.field}{Logical indicating whether or not a home field effect was modeled.}
\item{Hessian}{The Hessian of the model parameters, if requested.}
\item{parameters}{A vector of fitted model parameters.}

The function  \code{\link{game.pred}} may be used to predict the outcome of future games.
}
\references{

Karl, A.T., Yang, Y. and Lohr, S. (2013). Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to 
    Educational Value-Added Assessments. \emph{Computational Statistics and Data Analysis},  59, 13-27.

Karl, A., Yang, Y. and Lohr, S. (2014) Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects. \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Karl, A.T. (2012). The Sensitivity of College Football Rankings to Several Modeling Choices, \emph{Journal of Quantitative Analysis in Sports}, 
    Volume 8, Issue 3, DOI 10.1515/1559-0410.1471
}
\author{
Andrew T. Karl \email{akarl@asu.edu}, Jennifer Broatch
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{game.pred}}
}
\examples{
data(nfl2012)
mvglmmRank(nfl2012,method="PB0",first.order=TRUE,verbose=TRUE,max.iter.EM=1)
\donttest{
result <- mvglmmRank(nfl2012,method="PB0",first.order=TRUE,verbose=TRUE)
print(result)
game.pred(result,home="Denver Broncos",away="Green Bay Packers")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

