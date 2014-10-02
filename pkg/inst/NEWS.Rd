\name{NEWS}
\title{News for Package 'stabs'}

\section{Changes in stabs version 0.1-0 (2014-10-02)}{
  \itemize{
    \item Initial submission to CRAN
    \item \pkg{stabs} implements stability selection for a range of
    models, including a novel \code{matrix} interface for \code{stabsel}
    that can be used with generic fit functions. 
    \item The code is based on the \code{stabsel} function that was
    implemented in \pkg{mboost} until version 2.3-0. \code{stabsel} is
    now a generic function with a specialized method for boosting models
    that is implemented in \pkg{mboost}. The back-end is now completely
    implemented in package \pkg{stabs}.
  }
}