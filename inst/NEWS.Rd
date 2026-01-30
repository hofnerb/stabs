\name{NEWS}
\title{News for Package 'stabs'}

\section{Changes in stabs version 0.7-1 (2026-01-30)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item{Fixed timing issues for checks of examples on CRAN: 
      added \code{\\donttest} in two cases and make sure that code is not
      executed in parallel.}
    }
  }
}

\section{Changes in stabs version 0.7-0 (2026-01-27)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Allow \code{cutoff} to be less than 0.5 in case of 
        \code{assumption = "r-concave"}
        (closes \href{https://github.com/hofnerb/stabs/issues/35}{#35}).
      \item Updated error messages and error handling if the required package 
        for \code{fitfun} is missing (closes \href{https://github.com/hofnerb/stabs/issues/29}{#29}).
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item{\code{glmnet.lasso_maxCoef} now correctly extracts the vatriables 
        with largest (absolute) coefficients (closes \href{https://github.com/hofnerb/stabs/issues/32}{#32}).}
      \item{Fixed issues from CRAN checks: changed citation style, 
      fixed call to \code{\\code} in manuals and removed \code{LazyData} from \file{DESCRIPTION}. 
      Additionally added a link to the github issue tracker under \code{BugReports}
      (contributed by Roberto Villegas-Diaz; see pull request \href{https://github.com/hofnerb/stabs/pull/34}{#34}).}
    }
  }
}

\section{Changes in stabs version 0.6-4 (2021-01-28)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Removed support for graphical models via package \CRANpkg{QUIC}
        as QUIC is no longer maintained and was archived on CRAN.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item{Fixed CRAN checks regarding \code{all.equal} by adding 
        \code{check.environment = FALSE}.}
    }
  }
}

\section{Changes in stabs version 0.6-3 (2017-07-19)}{
  \subsection{Miscellaneous}{
    \itemize{
    \item Added reference to stability selection with \CRANpkg{gamboostLSS}.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
    \item Make vignettes conditional on suggested packages 
      (closes \href{https://github.com/hofnerb/stabs/issues/24}{#24}).
    \item Fixed tests.
    }
  }
}

\section{Changes in stabs version 0.6-2 (2017-01-30)}{
  \subsection{Miscellaneous}{
    \itemize{
    \item Move packages from \code{Enhances} to \code{Suggests}.
    \item Make manuals work without package \CRANpkg{TH.data}.
    \item Speed up examples by adding \code{\\donttest} in two occassions.
    \item Use markup \code{\\doi} for DOIs in manual.
    }
  }
}

\section{Changes in stabs version 0.6-1 (2017-01-28)}{
  \subsection{Miscellaneous}{
    \itemize{
    \item Added DOIs to \file{DESCRIPTION} to fulfill CRAN policies.
    \item Updated \file{README.md} and vignettes.
    }
  }
}

\section{Changes in stabs version 0.6-0 (2017-01-27)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Added support for graphical models via package \CRANpkg{QUIC}
        (closes \href{https://github.com/hofnerb/stabs/issues/12}{#12}, 
        \href{https://github.com/hofnerb/stabs/issues/13}{#13} and 
        \href{https://github.com/hofnerb/stabs/issues/14}{#14}; Code contributed
        by Richard Beare).
      \item Better handling of errors in (single) folds of \code{stabsel}:
        results of folds without errors are used and a \code{warning} is
        issued  (closes \href{https://github.com/hofnerb/stabs/issues/2}{#2}
        and \href{https://github.com/hofnerb/stabs/issues/7}{#7}). 
      \item Parallel computing via \code{mclapply}: Set
        \code{mc.preschedule = FALSE} per default (closes
        \href{https://github.com/hofnerb/stabs/issues/1}{#1}).
      \item Make sure that per default \code{glmnet.lasso} selects at maximum
        \code{q} variables. The previous anti-conservative version 
        stays available via \code{args.fitfun = list(type = "anticonservative")} (closes 
        \href{https://github.com/hofnerb/stabs/issues/23}{#23})
      }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Updated references: Hofner et al (2015) now available at BMC
        Bioinformatics. See \code{citation("stabs")} for details.
      \item Updated \file{DESCRIPTION} to be more informative.
      \item Make \code{fitfun}s more quiet (closes 
        \href{https://github.com/hofnerb/stabs/issues/19}{#19};
        Code contributed by Gokcen Eraslan)
      \item Replaced \code{base-learner} with \code{effect} in warnings and 
        error messages (closes 
        \href{https://github.com/hofnerb/stabs/issues/10}{#10}).
      \item Added example on using \code{glmnet.lasso} (closes 
        \href{https://github.com/hofnerb/stabs/issues/4}{#4}).
      \item Boosting specific changes: warn if \code{mstop} is to small.  
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Disallow specification of penalty parameters via \code{args.fitfun}
        (closes \href{https://github.com/hofnerb/stabs/issues/22}{#22}).
      \item Fixed checks for results of \code{fitfun}s 
        (closes \href{https://github.com/hofnerb/stabs/issues/17}{#17}; Code contributed
        by Andrey Tovchigrechko)
      \item Fixed citation. (closes \href{https://github.com/hofnerb/stabs/issues/9}{#9}).
      \item Fix plot labels if matrices are used (as opposed to \code{data.frame}s)
        (closes \href{https://github.com/hofnerb/stabs/issues/6}{#6}).
      \item Fix issue when variables are dropped from active set in 
        \code{lars.lasso} (closes \href{https://github.com/hofnerb/stabs/issues/5}{#5}).
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
      \CRANpkg{mboost}); \code{selected()} can now be used on stability
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
      \item Package development moved from R-forge
      to \url{https://github.com/hofnerb/stabs}

      \item Added \code{specifiedPFER} to \code{stabsel} and
      \code{stabsel_parameter} objects (needed for better outputs)
      \item Added \file{inst/CITATION}
      \item Added reference to stability selection paper
      \item \CRANpkg{stabs} now enhances \CRANpkg{gamboostLSS}
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
    \item \CRANpkg{stabs} implements stability selection for a range of
    models, including a novel \code{matrix} interface for \code{stabsel}
    that can be used with generic fit functions.
    \item The code is based on the \code{stabsel} function that was
    implemented in \CRANpkg{mboost} until version 2.3-0. \code{stabsel} is
    now a generic function with a specialized method for boosting models
    that is implemented in \CRANpkg{mboost}. The back-end is now completely
    implemented in package \CRANpkg{stabs}.
  }
}