\name{f2011}
\alias{f2011}
\docType{data}
\title{
2011 FBS College Football Regular Season Data
}
\description{
2011 FBS College Football Regular Season Data
}
\usage{data(f2011)}
\format{
  A data frame with 781 observations on the following 7 variables.
  \describe{
    \item{\code{home}}{a factor}
    \item{\code{Game.Date}}{a POSIXlt}
    \item{\code{away}}{a factor}
    \item{\code{home.score}}{a numeric vector}
    \item{\code{away.score}}{a numeric vector}
    \item{\code{neutral.site}}{a numeric vector}
    \item{\code{partition}}{a numeric vector}
  }
}
\source{
http://web1.ncaa.org/mfb/download.jsp?year=2011&div=IA
}
\examples{
data(f2011)
## maybe str(f2011) ; plot(f2011) ...
}
\keyword{datasets}
