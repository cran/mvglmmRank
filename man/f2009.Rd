\name{f2009}
\alias{f2009}
\docType{data}
\title{
2009 FBS College Football Regular Season Data
}
\description{
2009 FBS College Football Regular Season Data
}
\usage{data(f2009)}
\format{
  A data frame with 772 observations on the following 7 variables.
  \describe{
    \item{\code{home}}{a factor}
    \item{\code{Game.Date}}{a POSIXlt date variable}
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
http://web1.ncaa.org/mfb/download.jsp?year=2009&div=IA
}
\examples{
data(f2009)
## maybe str(f2009) ; plot(f2009) ...
}
\keyword{datasets}
