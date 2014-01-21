\name{export.nodal}
\alias{export.nodal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Outputs all node-level measures computed by the user
}
\description{
    Takes any output from nodal functions (\code{niches}, \code{nodal.local}, \code{nodal.global}, or \code{nodal.network}) and presents it in one matrix for further analysis outside of Blaunet or output to disk.
}
\usage{
export.nodal(blauObj, niches = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blauObj}{
    An object of class blau initialized with the function blau. User must have called at least one of the following functions: \code{niches}, \code{nodal.local}, \code{nodal.global}, or \code{nodal.network}. Otherwise, there is nothing to export.
}
  \item{niches}{
    Defaults to \code{TRUE}. If \code{TRUE}, includes the full matrix of individuals located in niches as part of the output. This information may not always be desirable for output, and setting this to \code{FALSE} will prevent the full matrix of individuals located in niches from being included.
}
}
\details{
    This function is useful to manipulate the computed nodal Blau statuses and analyze them outside the blaunet program, either within R or by exporting them to another statistical package. 
}
\value{
    Returns an object of class data.frame that includes all output from \code{niches}, \code{nodal.local}, \code{nodal.global}, or \code{nodal.network} previously computed by the user. The user must manually call one or more of these four functions prior to calling \code{export.nodal}. Columns are labeled with the name of the measure. Row names are node names provided with the blau function.

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
#blau object will store whatever you compute
b <- nodal.global(b)
b <- nodal.local(b, focal.niche = 'grppta')
export.nodal(b) #will export global and local
export.nodal(b, niches = FALSE) #suppress niche export

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
