\name{ncaab2012}
\alias{ncaab2012}
\docType{data}
\title{
2012 NCAA Division I Basketball Results
}
\description{
2012 NCAA Division I Basketball Results
}
\usage{data(ncaab2012)}
\format{
  A data frame with 5253 observations on the following 8 variables.
  \describe{
    \item{\code{date}}{a factor}
    \item{\code{away}}{a factor}
    \item{\code{home}}{a factor}
    \item{\code{neutral.site}}{a numeric vector}
    \item{\code{partition}}{a numeric vector}
    \item{\code{home_win}}{a numeric vector}
    \item{\code{home.score}}{a numeric vector}
    \item{\code{away.score}}{a numeric vector}
  }
}
\source{
http://masseyratings.com/data.php
}
\examples{
data(ncaab2012)
## maybe str(ncaab2012) ; plot(ncaab2012) ...
}
\keyword{datasets}
