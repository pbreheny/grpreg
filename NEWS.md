# grpreg 3.3.1 (2021-03-26)
  * Fixed: AUC() now compatible with survival 3.2.10
  * Internal: Fixed memory leak

# grpreg 3.3.0 (2020-06-10)
  * Fixed: sqrt(K) no longer hard-coded into discarding rules (thank you to Dan
    Kessler for pointing this out)
  * Testing: Now uses the tinytest package
  * Documentation: Removing references to grpregOverlap (hope to merge)

# grpreg 3.2.2 (2020-02-14)
  * Change: Better error detection for ill-conditioned, unpenalized matrices
  * Fixed: loss.grpsurv now works for total=FALSE
  * Internal: Lots of internal changes for cleaner, more reliable code
  * New version numbering system

# grpreg 3.2-1 (2019-02-26)
  * Change: Cross-validation now balances censoring across folds for survival
    models
  * Fixed: Leave-one-out cross-validation now works correctly for logistic
    regression

# grpreg 3.2-0 (2018-09-27)
  * New: cv.grpsurv now calculates SE, with bootstrap option
  * Change: R^2 now consistently uses the Cox-Snell definition for all types
    of models
  * Change: Survival loss now uses deviance
  * Change: cv.grpsurv now uses 'fold', not 'cv.ind', to declare assignments
  * Fixed: cv.grpreg now correctly handles out-of-order groups for Poisson
  * Fixed: cv.grpsurv now correctly standardizes out-of-order groups
  * Fixed: grpreg no longer returns loss=NA with family='binomial' for some
    lambda values
  * Internal: SSR-BEDPP optimization reinstated after bug fix
  * Internal: C code for binom/pois combined into gdfit_glm, lcdfit_glm
  * Documentation: Lots of updates
  * Documentation: vignette now html (used to be pdf)
  * Documentation: pkgdown website

# grpreg 3.1-4 (2018-06-15)
  * Fixed: Works with arbitrarily "messy" group structures now (constant
    columns, out of order groups, etc.) due to restructuring of standardization/
    orthogonalization
  * Internal: SSR-BEDPP rule turned off due to bug

# grpreg 3.1-3 (2018-04-07)
  * Internal: C code now uses || instead of |

# grpreg 3.1-2 (2017-07-05)
  * Fixed: Bug in applying screening rules with group lasso for linear
    regression with user-specified lambda sequence (thank you very much to
    Natasha Sahr for pointing this out)

# grpreg 3.1-1 (2017-06-07)
  * Fixed: Cross-validation no longer fails when constant columns are present
    (thank you to Matthew Rosenberg for pointing this out)
  * Fixed: Cross-validation no longer fails when group.multiplier is specified

# grpreg 3.1-0 (2017-05-18)
  * New: Additional tests and support for coersion of various types with
    respect to both X and y
  * Change: Convergence criterion now based on RMSD of linear predictors
  * Change: 'Lung' and 'Birthwt' data sets now use factor representation of
    group, as character vectors are inherently ambiguous with respect to order
  * Change: max.iter now based on total number of iterations for entire path
  * Internal: 'X', 'group', and 'group.multiplier' now bundled together in an
    object called 'XG' to enforce agreement at all times
  * Internal: new SSR-BEDPP feature screening rule for group lasso
  * Internal: Registration of native routines
  * Internal: Changing PROTECT/UNPROTECT to conform to new coding standards
  * Fixed: The binding of X and G fixes several potential bugs, including
    Issue #12 (GitHub)

# grpreg 3.0-2
  *  Fixed bug involving mismatch between group.multiplier and group if
	   group is given out of order.

# grpreg 3.0-1 (2016-06-06)
  * Fixed: memory allocation bug
  * Deprecation: Re-introduced 'birthwt.grpreg' for backwards compatibility,
    but this is deprecated

# grpreg 3.0-0 (2016-06-02)
  * New: methods for survival analysis (Cox modeling): grpsurv, cv.grpsurv, AUC,
    predict.grpsurv
  * New: option to return fitted values from cross-validation folds
    (returnY=TRUE) to cv.grpreg and cv.grpsurv
  * New: Added user interrupt checking
  * Change: Reformatted (and renamed) example data set 'Birthwt'; added example
    data set 'Lung' for survival
  * Internal: Greatly expanded suite of tests; various bugs identified and fixed
    as a result
  * Documentation: Added vignettes (a quick-start guide and a detailed
    description of available penalties)

# grpreg 2.8-1 (2015-05-30)
  * New: cv.grpreg now allows user to specify lambda (thanks to Vincent
    Arel-Bundock for suggesting this change)
  * Fixed: bug for predict.grpreg(fit, type="nvars") or type="ngroups" when
    scalar lambda value is passed
  * Documentation: Updated citations

# grpreg 2.8-0 (2014-11-15)
  * New: More flexible interface through the 'group' argument; groups may now be
    out of order, and may be named rather than only consecutive integers
  * New: 'X' can now be a matrix of integers (previously this would result in
    the passing of an incompatible storage type to C)
  * New: Additional error checks to prevent cryptic error messages
  * Internal: modifications to convergence monitoring
  * New: Added corrected AIC and extended BIC as options with select()
  * Change: summary.cv.grpreg now describes multitask learning models more
    accurately
  * Fixed: bug for multitask learning when number of outcomes = 2 (thank you to
    Aluma Dembo for pointing this out)
  * Fixed: Cross-validation for multitask learning now respects the multivariate
    structure of the response matrix
  * Fixed: bug in cv.grpreg when attempting to use leave-one-out
    cross-validation

# grpreg 2.7-1 (2014-08-13)
  * Fixed: More rigorous initialization at C level to prevent possible memory
    access problems
  * Fixed: predict() for types 'vars', 'nvars', and 'ngroups' with multivariate
    outcomes
  * Fixed: As a consequence of the above fix, summary(cvfit) now works for
    multivariate outcomes (thank you to Cajo ter Braak for pointing out that
    this was broken)

# grpreg 2.7-0 (2014-08-13)
  * New: support for Poisson regression
  * Internal: .Call now used instead of .C
  * Fixed: bug in cv.grpreg when attempting to use leave-one-out
    cross-validation (thank you to Cajo ter Braak for pointing this out)

# grpreg 2.6-0 (2014-03-21)
  * Internal: Various internal changes to make the package more efficient for
    large data sets

# grpreg 2.5-0 (2013-12-24)
  * New: group exponential lasso 'gel' method
  * New: 'gmax' option
  * New: 'nvars' and 'ngroups' options for predict
  * Change: appearance of summary.cv.grpreg display

# grpreg 2.4-0 (2013-06-07)
  * New: options in plot.cv.grpreg to plot estimates of r-squared,
    signal-to-noise ratio, scale parameter, and prediction error in addition to
    cross-validation error (deviance)
  * New: grpreg and cv.grpreg now allow matrix y to facilitation group penalized
    methods for seemingly unrelated regressions/multitask learning.  This is
    something of a 'beta' release at this point, and will be developed and
    refined further in future releases.
  * New: 'summary' method for cv.grpreg objects
  * New: 'coef' and 'predict' methods for cv.grpreg objects
  * Change: Brought gBridge up to date so that it now handles constant columns,
    etc. (see # grpreg 2.2-0)
  * Fixed: bug in predict type='coefficients' when 'lambda' argument specified
  * Fixed: bug in cv.grpreg with user-defined lambda values

# grpreg 2.3-0 (2013-02-10)
  * Internal: Switched to SVD-based orthogonalization to allow for linear
    dependency within groups

# grpreg 2.2-1 (2012-11-16)
  * Fixed: compilation error for 32-bit Windows
  * Fixed: bug in calculation of binomial deviance when fitted probabilities
    are close to 0 or 1

# grpreg 2.2-0 (2012-10-09)
  * New: select now Now allows '...' options to be passed to logLik
  * New: Added option to plot norm of each group, rather than individual
    coefficients
  * New: 'vars', 'groups', and 'norm' options added to 'predict'
  * Change: cv.grpreg now returns full data fit as well as CV errors; this
    allows cv.grpreg to handle constant columns and fixes some bugs
  * Fixed: logLik no longer calculates (meaningless) log-likelihoods for
    saturated models (thank you to Xiaowei Ren for pointing this out)
  * Fixed: bug for returning group when some groups were eliminated due to
    constant columns

# grpreg 2.1-0 (2012-07-28)
  * New: grpreg can now handle constant columns (they produce beta=0)
  * Fixed: Bug involving orthogonalization with unpenalized groups
  * Internal: restructuring of C code

# grpreg 2.0-0 (2012-07-21)
  * New: Group MCP, group SCAD methods added
  * New: Added 'cv.grpreg' to facilitate cross-validation
  * New: 'dfmax' option
  * New: 'group.multiplier' option
  * New: Allows specification of unpenalized groups
  * Change: gBridge now divorced from grpreg and given separate function
  * Internal: New algorithm for group lasso
  * Internal: Extensive internal refactoring of code
  * Internal: standardize and orthogonalize functions added
  * Internal: Much more extensive and reproducible code testing

# grpreg 1.2-0 (2011-06-22)
  * New: grpreg now returns 'loss'
  * New: Added logLik method
  * Change: Syntax of 'select' modified (no longer requires X, y to be passed)
  * Change: 'plot.grpreg' function more flexible
  * Change: 'n.lambda' to 'nlambda' in grpreg
  * Change: 'a' to 'gamma' for MCP tuning parameter
  * Change: 'lambda2' to 'alpha'
  * Removed: 'monitor' no longer an option in grpreg
  * Removed: 'criteria' option for select
  * Fixed: Bug in calculation of df for gLasso (grpreg.c)
  * Documendation: Updated citation and contact information
