\name{NEWS}
\title{News for Package 'stabs'}

\section{Changes in stabs version 0.5-2 (2015-0x-xx)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Better handling of errors in (single) folds of \code{stabsel}:
      results of folds without errors are used and a \code{warning} is
      issued  (closes
      \href{https://github.com/hofnerb/stabs/issues/2}{#2}). 
      \item Parallel computing via \code{mclapply}: Set
      \code{mc.preschedule = FALSE} per default (closes
      \href{https://github.com/hofnerb/stabs/issues/1}{#1}).
      \item Boosting specific changes: warn if \code{mstop} is to small.
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Updated references: Hofner et al (2015) now available at BMC
      Bioinformatics. See \code{citations("stabs")} for details.
      \item Updated \file{DESCRIPTION} to be more informative.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Adhere to CRAN policies regarding import of base packages
      (closes \href{https://github.com/hofnerb/stabs/issues/3}{#3}). 
    }
  } 
}

\section{Changes in stabs version 0.5-1 (2015-01-12)}{
  \subsection{Bug-fixes}{
    \itemize{
      \item Changes in \file{inst/CITATION} to make CRAN happy:
      Citations can now be extracted without the need to install the
      package.
      \item Changed \code{require} to \code{requireNamespace}
    }
  }
}

\section{Changes in stabs version 0.5-0 (2014-12-09)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Added a function \code{stabsel.stabsel()} to compute (new)
      parameter combinations for a fitted \code{stabsel} object
      \item Included \code{selected()} method (originally from
      \pkg{mboost}); \code{selected()} can now be used on stability
      selection results to extract selected effects
      \item Added functionality to extract parameters from fitted
      \code{stabsel} objects:
      use \code{parameters()} or \code{stabsel_parameters()}
      \item Export \code{run_stabsel} to make the actual \dQuote{fit}
      function accessible for other packages.
      \item Tweaked \code{print} methods a bit
      \item Added output that states the significance level.
      \item Added new lasso based fit (\code{glmnet.lasso_maxCoef})
      function where the \code{q} strongest predictors (according to the
      coefficient size) are selected (feature request from Rajen Shah
      <R.Shah _at_ statslab.cam.ac.uk>)
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Package development moved from
      \url{http://r-forge.r-project.org/projects/stabsel/}
      to \url{http://github.com/hofnerb/stabs}

      \item Added \code{specifiedPFER} to \code{stabsel} and
      \code{stabsel_parameter} objects (needed for better outputs)
      \item Added \file{inst/CITATION}
      \item Added reference to stability selection paper
      \item \pkg{stabs} now enhances \pkg{gamboostLSS}
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Made tests conditional on availability of packages
    }
  }
}

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