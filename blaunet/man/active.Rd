\name{active}
\alias{active}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Quick summary of blau object.
}
\description{
    Identifies measures that have been computed.
}
\usage{
active(blauObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blauObj}{
    An object of class blau initialized with the function \code{blau}.}
}
\details{
    Since measures computed with many functions are stored in the blau object itself, it is helpful to quickly see which elements are active.
}
\value{
A string indicating the active elements.

}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
data(TwoCities)
b <- blau(TwoCities, node.ids = 'respID', ecology.ids = 'samp')
b <- niches(b)
active(b) #will produce command line output

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
