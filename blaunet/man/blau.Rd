\name{blau}
\alias{blau}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts raw data into an object for Blau status analysis
}
\description{
Converts a matrix of organization memberships and demographic dimensions, along with (optionally) an edgelist or adjacency matrix, into an object of class \code{blau}. Automatically detects organization membership and demographic columns.}
\usage{
blau(
  square.data, graph = NULL, directed.el = FALSE, node.ids = NULL,
  weights = NULL, ecology.ids = NULL, exclude = NULL, dimensions = NULL,
  memberships = NULL, complete.cases = FALSE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{square.data}{
    R matrix or \code{data.frame} object that must contain demographic and membership information (in columns). May also include columns of individual or ecology identifiers, weights, or a primary membership column.}
  \item{graph}{
    A named edgelist or adjacency matrix. This is required for computing measures which incorporate network information. Relies on the \code{network} object from package \code{network}.
}
  \item{directed.el}{
    Defaults to FALSE. Used only to indicate whether an edgelist passed to this function should be treated as directed or undirected. Not necessary if passing an adjacency matrix.
}
  \item{node.ids}{
    Indicates the column which holds node (individual) identifiers. May be the column number or column name.
}
  \item{weights}{
    A column with weights for each observation. May be the column number or column name. In the absence of specification, all weights are assumed to be equal (and are set to 1). Weights are used in a weighted standard deviation calculation.
}
  \item{ecology.ids}{
    Indicates the column which holds ecology identifiers. May be the column number or column name.
}
  \item{exclude}{
    A way to manually exclude columns from automatic incorporation as membership or demographic columns. May be a vector of column names as strings, or a vector of column numbers. Useful for larger datasets where the vast majority of columns are included.
}
  \item{dimensions}{
    Indicates the columns which hold Blau parameter information. May be a vector of column names as strings, or a vector of column numbers. In the absence of specification, all non-binary columns that are not used for other purposes will be assumed to be demographic variables.
}
  \item{memberships}{
    Indicates the columns which hold organizational membership information.  May be a vector of column names as strings, or a vector of column numbers. In the absence of specification, all binary columns that are not used for other purposes will be assumed to be membership indicators.
}
  \item{complete.cases}{
    Defaults to \code{FALSE}. A boolean setting indicating whether all rows with at least one missing value should be dropped before proceeding with calculations. Otherwise, the program will be ``greedy'' about determining niche boundaries, using as much information as possible.}
}
\details{
The most common raw data format is a \code{.csv} file that contains ecology identifiers, node identifiers, blau parameters, and memberships, among other variables

The vast majority of configuration takes place when calling the blau function. As such, it is essential that the user understand how choices made here impact the operation of subsequent functions. The easiest way to get started is to determine which of the four optional parameters---node identifiers (\code{node.ids}), ecology identifiers (\code{ecology.ids}), weights (\code{weights}), and relational data (\code{graph})---are present in your dataset and will be used for analysis. These should be specified by indicating their locations with the appropriate function argument, and the \code{blau} function will automatically assume all other columns are either membership or demographic columns. If there are columns to be excluded from analysis, they can be specified with the exclude parameter. This type of setup is appropriate for the vast majority of datasets.

It is important to remember that any information incorporated into the \code{blau} object through this function will be used when calling subsequent functions. For instance, if your analysis does not require weights, but they are provided in the dataset, they should be explicitly excluded with the \code{exclude} argument.

If ecology identifiers are provided, all subsequent analyses will automatically proceed on a by-ecology level (unless specified explicitly in subsequent functions).

With network information, the most important consideration is that node identifiers are properly indicated and may be matched up with node identifications provided with the \code{node.ids} parameter. Adjacency matrix or edgelist input formats are both converted to an \code{network} object. The preferred format is a named edgelist (two columns, with node names in each row indicating an edge). 

If node names are \code{numeric} rather than \code{character}, they should still be specified in the input function with \code{node.ids} and a network should correctly indicate node identifiers.

If complete.cases is \code{FALSE} (the default option), we automatically use as much information as possible to compute niche boundaries. For example, an individual may have missing Blau parameter information for a certain dimension. Under the default settings, we still utilize the user's other demographic information to compute niche boundaries. If compleCases is specified as \code{TRUE}, then only observations with no missing values along all elements in the input matrix will be utilized in determining boundaries.
}
\value{
    \item{obj}{An object of class \code{blau}. Object contains several elements that may be accessed with the $ operator.}
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
##simple example
data(TwoCities)
b <- blau(TwoCities, node.ids = 'respID', ecology.ids = 'samp')
##blau object will store whatever you compute
b <- niches(b)
##see active elements
print(b)
##compute global Blau statuses
b <- nodal.global(b)

##assume we don't care about the 'ideo' column
b <- blau(TwoCities, node.ids = 'respID', ecology.ids = 'samp', exclude='ideo')
##compute niches like before
b <- niches(b)

##example with relational data
data(BSANet)
square.data <- BSANet$square.data
el <- BSANet$el

b <- blau(square.data, node.ids = 'person', ecology.ids = 'city', graph = el)
#compute dyadic statuses
b <- dyadic(b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
