\name{dyadic}
\alias{dyadic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Computes dyadic Blau status measures.
}
\description{
    Computes six dyadic measures (computed for each edge in the edgelist). They are Euclidean Distance, Mahalanobis Distance, Co-nichers, Co-outsiders, Straddlers, Spanners.
}
\usage{
dyadic(blauObj, dev.range = 1.5, ecologies.off = FALSE, m.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blauObj}{
    An object of class blau initialized with the function \code{blau}. Individuals will automatically be placed in niches with the function \code{niches} if this has not been done manually.
}
  \item{dev.range}{
    Defaults to 1.5. When creating niches, indicates standard deviation around the mean in each dimension to include in niche. A larger value will make niches larger and therefore include more individuals.}
  \item{ecologies.off}{
    Defaults to \code{FALSE}. If set to \code{TRUE}, treats all individuals as in the same ecology, even if ecology identifiers (\code{ecology.ids}) have been specified. Will call the niches function and overwrite its output even if it has been manually called by the user.
}
  \item{m.dist}{
    Defaults to \code{FALSE}. If \code{TRUE}, Mahalanobis distance will be computed in addition to Euclidean distance.
    }
}
\details{
All measures are at the edge-level of analysis and characterize ties/dyads rather than individuals/nodes.

Euclidean Distance: is the euclidean distance between two nodes in the edge

Mahalanobis Distance: is euclidean distance standardized by the unit of measurement

Co-nichers: are ties that are located within the same niche. A 0 indicates that both nodes of the tie do not share any niche in common. A 1 indicates that both nodes of the tie belong to only one common niche. A number of 2 or more indicates that the nodes of the tie share more than 1 common niche. 

Co-outsiders: is an indicator variable that denotes ties where both nodes do not belong to any niche. 

Straddlers: is an indicator variable that denotes ties where one node belongs to at least one niche and the other node does not belong to any niche. 

Spanners: is an indicator variable that denotes ties where one node belongs to one niche and the other node belongs to a different niche.
}
\value{
    Returns a matrix stored in object$dyadic that contains eight columns. The first two contain the ego and alter for the edge, and the next six contain the dyadic measures.

}
%%}
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
el <- BSANet$el #edgelist
adj <- BSANet$adj #adjacency matrix

b <- blau(square.data, node.ids = 'person', ecology.ids = 'city', graph = el)
#implicitly a directed edgelist
#will automatically compute niches
b <- dyadic(b)

#adjacency matricies are also OK
b <- blau(square.data, node.ids = 'person', ecology.ids = 'city', graph = adj)
#will automatically compute niches
b <- dyadic(b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
