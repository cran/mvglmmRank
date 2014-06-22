\name{nfl2012}
\alias{nfl2012}
\docType{data}
\title{
2012 NFL Regular Season Data
}
\description{
2012 NFL Regular Season Data
}
\usage{data(nfl2012)}
\format{
  A data frame with 256 observations on the following 7 variables.
  \describe{
    \item{\code{Date}}{a factor}
    \item{\code{away}}{a factor} 
    \item{\code{home}}{a factor} 
    \item{\code{neutral.site}}{a numeric vector}
    \item{\code{home.score}}{a numeric vector}
    \item{\code{away.score}}{a numeric vector}
    \item{\code{partition}}{a numeric vector}
  }
}
\source{
http://masseyratings.com/data.php
}
\examples{
data(nfl2012)
## maybe str(nfl2012) ; plot(nfl2012) ...
}
\keyword{datasets}
