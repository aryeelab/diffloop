% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loopFunctions.R
\docType{methods}
\name{removeSelfLoops}
\alias{removeSelfLoops}
\alias{removeSelfLoops,loops-method}
\title{Remove self loops}
\usage{
removeSelfLoops(dlo)

\S4method{removeSelfLoops}{loops}(dlo)
}
\arguments{
\item{dlo}{A loops object}
}
\value{
A loops object
}
\description{
\code{removeSelfLoops} removes instances where a loop is observed
between the same anchor
}
\details{
This function removes loops from the \code{interactions} slot that reference
the same index of the \code{anchors} slot.
}
\examples{
rda<-paste(system.file('rda',package='diffloop'),'loops.small.rda',sep='/')
load(rda)
jpn_unique <- removeSelfLoops(loops.small)

}

