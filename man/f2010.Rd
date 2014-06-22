\name{f2010}
\alias{f2010}
\docType{data}
\title{
2010 FBS College Football Regular Season Data
}
\description{
2010 FBS College Football Regular Season Data
}
\usage{data(f2010)}
\format{
  A data frame with 770 observations on the following 7 variables.
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
http://web1.ncaa.org/mfb/download.jsp?year=2010&div=IA
}
\examples{
data(f2010)
## maybe str(f2010) ; plot(f2010) ...
}
\keyword{datasets}
