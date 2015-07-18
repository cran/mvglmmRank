\name{game.pred}
\alias{game.pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict outcomes of games.
}
\description{
After fitting a model with \code{\link{mvglmmRank}}, \code{game.pred} uses that model to predict outcomes of future matchups.
}
\usage{
game.pred(res, home, away, neutral.site = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
an object of class mvglmmRank
}
  \item{home}{
a character string for the home team (use quotation marks!)
}
  \item{away}{
a character string for the away team (use quotation marks!)
}
  \item{neutral.site}{
logical. If TRUE, uses the neutral site mean score, assuming some of the games in the training data occured at neutral sites.
}
}
\value{
Prints predicted scores and/or predicted probability of a home team win, depending on the type of model specified by \code{res}.
}
\references{
Karl, A.T., Broatch, J. (2014). mvglmmRank: An R Package Providing Multivariate Generalized Linear Mixed Models for Evaluating Sporting Teams. \emph{Submitted.} 

Karl, A.T., Yang, Y. and Lohr, S. (2013). Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to 
    Educational Value-Added Assessments. \emph{Computational Statistics and Data Analysis},  59, 13-27.

Karl, A., Yang, Y. and Lohr, S. (2014) Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Karl, A.T. (2012). The Sensitivity of College Football Rankings to Several Modeling Choices, \emph{Journal of Quantitative Analysis in Sports}, 
    Volume 8, Issue 3, DOI 10.1515/1559-0410.1471
}
\author{
Andrew T. Karl and Jennifer Broatch
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


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

