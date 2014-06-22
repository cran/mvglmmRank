\name{nba2013}
\alias{nba2013}
\docType{data}
\title{
2013 NBA Data
}
\description{
2013 NBA Data
}
\usage{data(nba2013)}
\format{
  A data frame with 1229 observations on the following 9 variables.
  \describe{
    \item{\code{Date}}{a factor}
    \item{\code{away}}{a factor}
    \item{\code{home}}{a factor}
    \item{\code{OT}}{a factor}
    \item{\code{partition}}{a numeric vector}
    \item{\code{neutral.site}}{a numeric vector}
    \item{\code{ot.count}}{a numeric vector}
    \item{\code{home.score}}{a numeric vector}
    \item{\code{away.score}}{a numeric vector}
  }
}
\source{
http://masseyratings.com/data.php
}
\examples{
data(nba2013)
## maybe str(nba2013) ; plot(nba2013) ...
}
\keyword{datasets}
