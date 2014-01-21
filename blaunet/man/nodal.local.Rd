\name{nodal.local}
\alias{nodal.local}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Computes Blau statuses for individuals with respect to a primary membership
}
\description{
    Computes two measures: whether an individual is within the niche of the primary membership specified, and whether the individual is a member of the primary membership but outside of that membership's niche.
}
\usage{
nodal.local(blauObj, focal.niche = NULL, dev.range = 1.5, ecologies.off = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blauObj}{
    An object of class blau initialized with the function blau. Individuals will automatically be placed in niches with the function \code{niches} if this has not been done manually.
}
  \item{focal.niche}{
    Specifies a focial niche by name of the organization.
}
  \item{dev.range}{
    Defaults to 1.5. When creating niches, indicates standard deviation around the mean in each dimension to include in niche. A larger value will make niches larger and therefore include more individuals.
}
  \item{ecologies.off}{
    Defaults to \code{FALSE}. If set to \code{TRUE}, treats all individuals as in the same ecology, even if ecology identifiers (\code{ecology.ids}) have been specified. Will call the niches function and overwrite its output even if it has been manually called by the user.
}
}
\details{
This provides information about the focal membership for each individual.

FocNicher: Provides information on how many niches a person is in including the focal niche (it may also be accessed manually by \code{object$isInNiche}). A value of 0 indicates a person is not in any niche. A value of 1 indicates the person is only in the focal niche. A value of more than 1 indicates the person is in more than the focal niche. ((it may be more intuitive if we use two columns: in focal niche and total # niches))

MemNotNiche: indicates whether an individual is a member of the focal membership but outside of the focal niche. This indicates that an individual is atypical compared to other group members.

}
\value{
    A matrix stored in \code{object$nodalLocal} with two columns, each containing one of the two measures. Row names are node names provided with the \code{blau} function.
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
#will automatically compute niches
b0 <- nodal.local(b, focal.niche = 'grppta')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
