\name{f2012}
\alias{f2012}
\docType{data}
\title{
2012 FBS College Football Regular Season Data
}
\description{
2012 FBS College Football Regular Season Data
}
\usage{data(f2012)}
\format{
  A data frame with 809 observations on the following 9 variables.
  \describe{
    \item{\code{home}}{a factor}
    \item{\code{Game.Date}}{a POSIXlt}
    \item{\code{away}}{a factor}
    \item{\code{home.response}}{a numeric vector}
    \item{\code{home.score}}{a numeric vector}
    \item{\code{away.response}}{a numeric vector}
    \item{\code{away.score}}{a numeric vector}
    \item{\code{neutral.site}}{a numeric vector}
    \item{\code{partition}}{a numeric vector}
  }
}
\source{
http://web1.ncaa.org/mfb/download.jsp?year=2012&div=IA
}
\examples{
data(f2012)
## maybe str(f2012) ; plot(f2012) ...
}
\keyword{datasets}
