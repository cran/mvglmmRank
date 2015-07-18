\name{f2008}
\alias{f2008}
\docType{data}
\title{
2008 FBS College Football Regular Season Data
}
\description{
2008 FBS College Football Regular Season Data
}
\usage{data(f2008)}
\format{
  A data frame with 772 observations on the following 9 variables.
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
http://web1.ncaa.org/mfb/download.jsp?year=2008&div=IA
}
\examples{
data(f2008)
## maybe str(f2008) ; plot(f2008) ...
}
\keyword{datasets}
