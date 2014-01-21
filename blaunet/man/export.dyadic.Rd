\name{export.dyadic}
\alias{export.dyadic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Outputs dyadic level measures computed by the user
}
\description{
Takes any output from dyadic blau functions (dyadic) and presents it in one matrix for further analysis outside of Blaunet or output to disk.
}
\usage{
export.dyadic(blauObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blauObj}{
%%     ~~Describe \code{blauObj} here~~
}
}
\details{
This function is useful to manipulate the computed dyadic Blau statuses and analyze them outside the blaunet program, either within R or by exporting them to another statistical package. 
}
\value{
    Returns an object of class data.frame that includes all output from function dyadic. The user must manually call the \code{dyadic} function prior to calling \code{export.dyadic}.
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
data(BSANet)
square.data <- BSANet$square.data
el <- BSANet$el
b <- blau(square.data, node.ids = 'person', ecology.ids = 'city', graph = el)
b <- dyadic(b)
export.dyadic(b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
