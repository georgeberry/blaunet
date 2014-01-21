\name{nodal.network}
\alias{nodal.network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Computes nodal spanners
}
\description{
    Computes whether an individual spans (has connections in) another niche and the number of other niches each individual spans to.
}
\usage{
nodal.network(blauObj, dev.range = 1.5, ecologies.off = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blauObj}{
    An object of class blau initialized with the function \code{blau}. Individuals will automatically be placed in niches with the function \code{niches} if this has not been done manually.}
  \item{dev.range}{
    Defaults to 1.5. When creating niches, indicates standard deviation around the mean in each dimension to include in niche. A larger value will make niches larger and therefore include more individuals.}
  \item{ecologies.off}{
    Defaults to \code{FALSE}. If set to \code{TRUE}, treats all individuals as in the same ecology, even if ecology identifiers (\code{ecology.ids}) have been specified. Will call the \code{niches} function and overwrite its output even if it has been manually called by the user.}
}
\details{
This function requires sociometric network data and identifies ties that bridge niches. The preferred format is a named edgelist, although an adjacency matrix with properly named rows and columns is also acceptable. Network information will be matched to organization and demographic information by these names.

For individuals \code{i} and \code{j} in niches \code{alpha} and \code{beta}, if individual \code{i} is in niche \code{alpha} and individual \code{j} is in niche \code{beta}, and \code{i} and \code{j} have an edge between them, then \code{i} spans to niche \code{beta} and \code{j} spans to niche \code{alpha}. However, if \code{i} is in both niche \code{alpha} and niche \code{beta}, \code{i} would not span to \code{beta} because \code{i} is already in niche \code{beta}.
}
\value{
    A matrix stored in \code{object$nodalNetwork} with two columns holding the two measures. Row names are node names provided with the \code{blau} function.

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
el <- BSANet$el
square.data <- BSANet$square.data
b <- blau(square.data, node.ids = 'person', ecology.ids = 'city', graph = el)
#will automatically compute niches
b <- nodal.network(b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
