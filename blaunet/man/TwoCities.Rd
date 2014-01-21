\name{TwoCities}
\alias{TwoCities}
\docType{data}
\title{
Example social capital dataset
}
\description{
The dataset contains 1008 individuals, their attributes on 4 demographic characteristics for constructing Blau parameters, and their memberships in 18 voluntary organizations, as well as six other assorted variables. The dataset contains a sample of individuals from two cities in the United States: Bismarck, North Dakota and Grand Rapids, Michigan.
}
\usage{data(TwoCities)}
\format{
  A data frame with 1008 observations on the following 30 variables.
  \describe{
    \item{\code{samp}}{a numeric vector}
    \item{\code{respID}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{income}}{a numeric vector}
    \item{\code{educ}}{a numeric vector}
    \item{\code{wrktime}}{a numeric vector}
    \item{\code{grprel}}{a numeric vector}
    \item{\code{grpsport}}{a numeric vector}
    \item{\code{grpyouth}}{a numeric vector}
    \item{\code{grppta}}{a numeric vector}
    \item{\code{grpvet}}{a numeric vector}
    \item{\code{grpnei}}{a numeric vector}
    \item{\code{grpeld}}{a numeric vector}
    \item{\code{grpsoc}}{a numeric vector}
    \item{\code{grplab}}{a numeric vector}
    \item{\code{grpprof}}{a numeric vector}
    \item{\code{grpfrat}}{a numeric vector}
    \item{\code{grpeth}}{a numeric vector}
    \item{\code{grppol}}{a numeric vector}
    \item{\code{grpart}}{a numeric vector}
    \item{\code{grphob}}{a numeric vector}
    \item{\code{grpself}}{a numeric vector}
    \item{\code{grpwww}}{a numeric vector}
    \item{\code{grpothr}}{a numeric vector}
    \item{\code{gender}}{a numeric vector}
    \item{\code{race}}{a numeric vector}
    \item{\code{ideo}}{a numeric vector}
    \item{\code{trust}}{a numeric vector}
    \item{\code{friends}}{a numeric vector}
    \item{\code{divrsity}}{a numeric vector}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
The TwoCities dataset is excerpted with permission from the Social Capital Benchmark Survey, which was collected by Professor Robert D. Putnam of the Saguaro Seminar Civic Engagement in America, a project of the John F. Kennedy School of Government at Harvard University and numerous community foundations nationwide, and made available through the Roper Center for Public Opinion Research. The full data can be accessed at the Roper Center for Public Opinion website at the following url: \url{http://www.ropercenter.uconn.edu/data_access/data/datasets/social_capital_community_survey_2006.htm}
}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(TwoCities)
}
\keyword{datasets}
