<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1" />
<meta name="generator" content="pdoc 0.9.2" />
<title>dropkick.logistic API documentation</title>
<meta name="description" content="Logistic regression model functions for dropkick" />
<link rel="preload stylesheet" as="style" href="https://cdnjs.cloudflare.com/ajax/libs/10up-sanitize.css/11.0.1/sanitize.min.css" integrity="sha256-PK9q560IAAa6WVRRh76LtCaI8pjTJ2z11v0miyNNjrs=" crossorigin>
<link rel="preload stylesheet" as="style" href="https://cdnjs.cloudflare.com/ajax/libs/10up-sanitize.css/11.0.1/typography.min.css" integrity="sha256-7l/o7C8jubJiy74VsKTidCy1yBkRtiUGbVkYBylBqUg=" crossorigin>
<link rel="stylesheet preload" as="style" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/10.1.1/styles/github.min.css" crossorigin>
<style>:root{--highlight-color:#fe9}.flex{display:flex !important}body{line-height:1.5em}#content{padding:20px}#sidebar{padding:30px;overflow:hidden}#sidebar > *:last-child{margin-bottom:2cm}.http-server-breadcrumbs{font-size:130%;margin:0 0 15px 0}#footer{font-size:.75em;padding:5px 30px;border-top:1px solid #ddd;text-align:right}#footer p{margin:0 0 0 1em;display:inline-block}#footer p:last-child{margin-right:30px}h1,h2,h3,h4,h5{font-weight:300}h1{font-size:2.5em;line-height:1.1em}h2{font-size:1.75em;margin:1em 0 .50em 0}h3{font-size:1.4em;margin:25px 0 10px 0}h4{margin:0;font-size:105%}h1:target,h2:target,h3:target,h4:target,h5:target,h6:target{background:var(--highlight-color);padding:.2em 0}a{color:#058;text-decoration:none;transition:color .3s ease-in-out}a:hover{color:#e82}.title code{font-weight:bold}h2[id^="header-"]{margin-top:2em}.ident{color:#900}pre code{background:#f8f8f8;font-size:.8em;line-height:1.4em}code{background:#f2f2f1;padding:1px 4px;overflow-wrap:break-word}h1 code{background:transparent}pre{background:#f8f8f8;border:0;border-top:1px solid #ccc;border-bottom:1px solid #ccc;margin:1em 0;padding:1ex}#http-server-module-list{display:flex;flex-flow:column}#http-server-module-list div{display:flex}#http-server-module-list dt{min-width:10%}#http-server-module-list p{margin-top:0}.toc ul,#index{list-style-type:none;margin:0;padding:0}#index code{background:transparent}#index h3{border-bottom:1px solid #ddd}#index ul{padding:0}#index h4{margin-top:.6em;font-weight:bold}@media (min-width:200ex){#index .two-column{column-count:2}}@media (min-width:300ex){#index .two-column{column-count:3}}dl{margin-bottom:2em}dl dl:last-child{margin-bottom:4em}dd{margin:0 0 1em 3em}#header-classes + dl > dd{margin-bottom:3em}dd dd{margin-left:2em}dd p{margin:10px 0}.name{background:#eee;font-weight:bold;font-size:.85em;padding:5px 10px;display:inline-block;min-width:40%}.name:hover{background:#e0e0e0}dt:target .name{background:var(--highlight-color)}.name > span:first-child{white-space:nowrap}.name.class > span:nth-child(2){margin-left:.4em}.inherited{color:#999;border-left:5px solid #eee;padding-left:1em}.inheritance em{font-style:normal;font-weight:bold}.desc h2{font-weight:400;font-size:1.25em}.desc h3{font-size:1em}.desc dt code{background:inherit}.source summary,.git-link-div{color:#666;text-align:right;font-weight:400;font-size:.8em;text-transform:uppercase}.source summary > *{white-space:nowrap;cursor:pointer}.git-link{color:inherit;margin-left:1em}.source pre{max-height:500px;overflow:auto;margin:0}.source pre code{font-size:12px;overflow:visible}.hlist{list-style:none}.hlist li{display:inline}.hlist li:after{content:',\2002'}.hlist li:last-child:after{content:none}.hlist .hlist{display:inline;padding-left:1em}img{max-width:100%}td{padding:0 .5em}.admonition{padding:.1em .5em;margin-bottom:1em}.admonition-title{font-weight:bold}.admonition.note,.admonition.info,.admonition.important{background:#aef}.admonition.todo,.admonition.versionadded,.admonition.tip,.admonition.hint{background:#dfd}.admonition.warning,.admonition.versionchanged,.admonition.deprecated{background:#fd4}.admonition.error,.admonition.danger,.admonition.caution{background:lightpink}</style>
<style media="screen and (min-width: 700px)">@media screen and (min-width:700px){#sidebar{width:30%;height:100vh;overflow:auto;position:sticky;top:0}#content{width:70%;max-width:100ch;padding:3em 4em;border-left:1px solid #ddd}pre code{font-size:1em}.item .name{font-size:1em}main{display:flex;flex-direction:row-reverse;justify-content:flex-end}.toc ul ul,#index ul{padding-left:1.5em}.toc > ul > li{margin-top:.5em}}</style>
<style media="print">@media print{#sidebar h1{page-break-before:always}.source{display:none}}@media print{*{background:transparent !important;color:#000 !important;box-shadow:none !important;text-shadow:none !important}a[href]:after{content:" (" attr(href) ")";font-size:90%}a[href][title]:after{content:none}abbr[title]:after{content:" (" attr(title) ")"}.ir a:after,a[href^="javascript:"]:after,a[href^="#"]:after{content:""}pre,blockquote{border:1px solid #999;page-break-inside:avoid}thead{display:table-header-group}tr,img{page-break-inside:avoid}img{max-width:100% !important}@page{margin:0.5cm}p,h2,h3{orphans:3;widows:3}h1,h2,h3,h4,h5,h6{page-break-after:avoid}}</style>
<script defer src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/10.1.1/highlight.min.js" integrity="sha256-Uv3H6lx7dJmRfRvH8TH6kJD1TSK1aFcwgx+mdg3epi8=" crossorigin></script>
<script>window.addEventListener('DOMContentLoaded', () => hljs.initHighlighting())</script>
</head>
<body>
<main>
<article id="content">
<header>
<h1 class="title">Module <code>dropkick.logistic</code></h1>
</header>
<section id="section-intro">
<p>Logistic regression model functions for dropkick</p>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">&#34;&#34;&#34;
Logistic regression model functions for dropkick
&#34;&#34;&#34;
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.special import expit
from scipy.sparse import issparse, csc_matrix
from sklearn.base import BaseEstimator
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import check_array, check_X_y
from sklearn.utils.multiclass import check_classification_targets
from .errors import _check_error_flag

from _glmnet import lognet, splognet, lsolns
from .util import (
    _fix_lambda_path,
    _check_user_lambda,
    _interpolate_model,
    _score_lambda_path,
)


class LogitNet(BaseEstimator):
    &#34;&#34;&#34;
    Logistic Regression with elastic net penalty.
    This is a wrapper for the glmnet function lognet.

    Parameters
    ----------
    alpha : float, default 1
        The alpha parameter, 0 &lt;= alpha &lt;= 1, 0 for ridge, 1 for lasso

    n_lambda : int, default 100
        Maximum number of lambda values to compute

    min_lambda_ratio : float, default 1e-4
        In combination with `n_lambda`, the ratio of the smallest and largest
        values of lambda computed.

    lambda_path : array, default None
        In place of supplying n_lambda, provide an array of specific values
        to compute. The specified values must be in decreasing order. When
        None, the path of lambda values will be determined automatically. A
        maximum of `n_lambda` values will be computed.

    standardize : bool, default True
        Standardize input features prior to fitting. The final coefficients
        will be on the scale of the original data regardless of the value
        of standardize.

    fit_intercept : bool, default True
        Include an intercept term in the model.

    lower_limits : array, (shape n_features,) default -infinity
        Array of lower limits for each coefficient, must be non-positive.
        Can be a single value (which is then replicated), else an array
        corresponding to the number of features.

    upper_limits : array, (shape n_features,) default +infinity
        Array of upper limits for each coefficient, must be positive.
        See lower_limits.

    cut_point : float, default 1
        The cut point to use for selecting lambda_best.
            arg_max lambda  cv_score(lambda) &gt;= cv_score(lambda_max) - cut_point * standard_error(lambda_max)

    n_splits : int, default 3
        Number of cross validation folds for computing performance metrics and
        determining `lambda_best_` and `lambda_max_`. If non-zero, must be
        at least 3.

    scoring : string, callable or None, default None
        Scoring method for model selection during cross validation. When None,
        defaults to classification score. Valid options are `accuracy`,
        `roc_auc`, `average_precision`, `precision`, `recall`. Alternatively,
        supply a function or callable object with the following signature
        ``scorer(estimator, X, y)``. Note, the scoring function affects the
        selection of `lambda_best_` and `lambda_max_`, fitting the same data
        with different scoring methods will result in the selection of
        different models.

    n_jobs : int, default 1
        Maximum number of threads for computing cross validation metrics.

    tol : float, default 1e-7
        Convergence tolerance.

    max_iter : int, default 100000
        Maximum passes over the data.

    random_state : number, default None
        Seed for the random number generator. The glmnet solver is not
        deterministic, this seed is used for determining the cv folds.

    max_features : int
        Optional maximum number of features with nonzero coefficients after
        regularization. If not set, defaults to X.shape[1] during fit
        Note, this will be ignored if the user specifies lambda_path.

    verbose : bool, default False
        When True some warnings and log messages are suppressed.

    Attributes
    ----------
    classes_ : array, shape(n_classes,)
        The distinct classes/labels found in y.

    n_lambda_ : int
        The number of lambda values found by glmnet. Note, this may be less
        than the number specified via n_lambda.

    lambda_path_ : array, shape (n_lambda_,)
        The values of lambda found by glmnet, in decreasing order.

    coef_path_ : array, shape (n_classes, n_features, n_lambda_)
        The set of coefficients for each value of lambda in lambda_path_.

    coef_ : array, shape (n_clases, n_features)
        The coefficients corresponding to lambda_best_.

    intercept_ : array, shape (n_classes,)
        The intercept corresponding to lambda_best_.

    intercept_path_ : array, shape (n_classes, n_lambda_)
        The set of intercepts for each value of lambda in lambda_path_.

    cv_mean_score_ : array, shape (n_lambda_,)
        The mean cv score for each value of lambda. This will be set by fit_cv.

    cv_standard_error_ : array, shape (n_lambda_,)
        The standard error of the mean cv score for each value of lambda, this
        will be set by fit_cv.

    lambda_max_ : float
        The value of lambda that gives the best performance in cross
        validation.

    lambda_best_ : float
        The largest value of lambda which is greater than lambda_max_ and
        performs within cut_point * standard error of lambda_max_.
    &#34;&#34;&#34;

    CV = StratifiedKFold

    def __init__(
        self,
        alpha=1,
        n_lambda=100,
        min_lambda_ratio=1e-4,
        lambda_path=None,
        standardize=True,
        fit_intercept=True,
        lower_limits=-np.inf,
        upper_limits=np.inf,
        cut_point=1.0,
        n_splits=3,
        scoring=None,
        n_jobs=1,
        tol=1e-7,
        max_iter=100000,
        random_state=None,
        max_features=None,
        verbose=False,
    ):

        self.alpha = alpha
        self.n_lambda = n_lambda
        self.min_lambda_ratio = min_lambda_ratio
        self.lambda_path = lambda_path
        self.standardize = standardize
        self.lower_limits = lower_limits
        self.upper_limits = upper_limits
        self.fit_intercept = fit_intercept
        self.cut_point = cut_point
        self.n_splits = n_splits
        self.scoring = scoring
        self.n_jobs = n_jobs
        self.tol = tol
        self.max_iter = max_iter
        self.random_state = random_state
        self.max_features = max_features
        self.verbose = verbose

    def fit(self, adata, y, n_hvgs, sample_weight=None, relative_penalties=None):
        &#34;&#34;&#34;
        Fit the model to training data. If n_splits &gt; 1 also run n-fold cross
        validation on all values in lambda_path.

        The model will be fit n+1 times. On the first pass, the lambda_path
        will be determined, on the remaining passes, the model performance for
        each value of lambda. After cross validation, the attribute
        `cv_mean_score_` will contain the mean score over all folds for each
        value of lambda, and `cv_standard_error_` will contain the standard
        error of `cv_mean_score_` for each value of lambda. The value of lambda
        which achieves the best performance in cross validation will be saved
        to `lambda_max_` additionally, the largest value of lambda s.t.:
            cv_score(l) &gt;= cv_score(lambda_max_) -\
                           cut_point * standard_error(lambda_max_)
        will be saved to `lambda_best_`.

        Parameters
        ----------
        adata : anndata.AnnData, shape (n_samples, n_features)
            Input features in .X

        y : array, shape (n_samples,)
            Target values

        n_hvgs : int, number of highly-variable genes to feature-select

        sample_weight : array, shape (n_samples,)
            Optional weight vector for observations

        relative_penalties: array, shape (n_features,)
            Optional relative weight vector for penalty.
            0 entries remove penalty.

        Returns
        -------
        self : object
            Returns self.
        &#34;&#34;&#34;
        # X is scaled counts for all cells in HVGs only for first run
        X = adata[:, adata.var.highly_variable].X.copy()

        X, y = check_X_y(X, y, accept_sparse=&#34;csr&#34;, ensure_min_samples=2)
        if sample_weight is None:
            sample_weight = np.ones(X.shape[0])
        else:
            sample_weight = np.asarray(sample_weight)

        if not np.isscalar(self.lower_limits):
            self.lower_limits = np.asarray(self.lower_limits)
            if len(self.lower_limits) != X.shape[1]:
                raise ValueError(&#34;lower_limits must equal number of features&#34;)

        if not np.isscalar(self.upper_limits):
            self.upper_limits = np.asarray(self.upper_limits)
            if len(self.upper_limits) != X.shape[1]:
                raise ValueError(&#34;upper_limits must equal number of features&#34;)

        if (
            any(self.lower_limits &gt; 0)
            if isinstance(self.lower_limits, np.ndarray)
            else self.lower_limits &gt; 0
        ):
            raise ValueError(&#34;lower_limits must be non-positive&#34;)

        if (
            any(self.upper_limits &lt; 0)
            if isinstance(self.upper_limits, np.ndarray)
            else self.upper_limits &lt; 0
        ):
            raise ValueError(&#34;upper_limits must be positive&#34;)

        if self.alpha &gt; 1 or self.alpha &lt; 0:
            raise ValueError(&#34;alpha must be between 0 and 1&#34;)

        # fit the model
        self._fit(X, y, sample_weight, relative_penalties)

        # score each model on the path of lambda values found by glmnet and
        # select the best scoring
        if self.n_splits &gt;= 3:
            self._cv = self.CV(
                n_splits=self.n_splits, shuffle=True, random_state=self.random_state
            )

            cv_scores, _hvgs = _score_lambda_path(
                self,
                adata,
                y,
                n_hvgs,
                sample_weight,
                relative_penalties,
                self.scoring,
                n_jobs=self.n_jobs,
                verbose=self.verbose,
            )

            self.cv_mean_score_ = np.atleast_1d(np.mean(cv_scores, axis=0))
            self.cv_standard_error_ = np.atleast_1d(stats.sem(cv_scores))

            self.lambda_max_inx_ = np.argmax(self.cv_mean_score_)
            self.lambda_max_ = self.lambda_path_[self.lambda_max_inx_]

            target_score = (
                self.cv_mean_score_[self.lambda_max_inx_]
                - self.cut_point * self.cv_standard_error_[self.lambda_max_inx_]
            )

            self.lambda_best_inx_ = np.argwhere(self.cv_mean_score_ &gt;= target_score)[0]
            self.lambda_best_ = self.lambda_path_[self.lambda_best_inx_]

            self.coef_ = self.coef_path_[..., self.lambda_best_inx_]
            self.coef_ = self.coef_.squeeze(axis=self.coef_.ndim - 1)
            self.intercept_ = self.intercept_path_[..., self.lambda_best_inx_].squeeze()
            if self.intercept_.shape == ():  # convert 0d array to scalar
                self.intercept_ = float(self.intercept_)

        return self

    def _fit(self, X, y, sample_weight=None, relative_penalties=None):
        if self.lambda_path is not None:
            n_lambda = len(self.lambda_path)
            min_lambda_ratio = 1.0
        else:
            n_lambda = self.n_lambda
            min_lambda_ratio = self.min_lambda_ratio

        check_classification_targets(y)
        self.classes_ = np.unique(y)  # the output of np.unique is sorted
        n_classes = len(self.classes_)
        if n_classes &lt; 2:
            raise ValueError(&#34;Training data need to contain at least 2 &#34; &#34;classes.&#34;)

        # glmnet requires the labels a one-hot-encoded array of
        # (n_samples, n_classes)
        if n_classes == 2:
            # Normally we use 1/0 for the positive and negative classes. Since
            # np.unique sorts the output, the negative class will be in the 0th
            # column. We want a model predicting the positive class, not the
            # negative class, so we flip the columns here (the != condition).
            #
            # Broadcast comparison of self.classes_ to all rows of y. See the
            # numpy rules on broadcasting for more info, essentially this
            # &#34;reshapes&#34; y to (n_samples, n_classes) and self.classes_ to
            # (n_samples, n_classes) and performs an element-wise comparison
            # resulting in _y with shape (n_samples, n_classes).
            _y = (y[:, None] != self.classes_).astype(np.float64, order=&#34;F&#34;)
        else:
            # multinomial case, glmnet uses the entire array so we can
            # keep the original order.
            _y = (y[:, None] == self.classes_).astype(np.float64, order=&#34;F&#34;)

        # use sample weights, making sure all weights are positive
        # this is inspired by the R wrapper for glmnet, in lognet.R
        if sample_weight is not None:
            weight_gt_0 = sample_weight &gt; 0
            sample_weight = sample_weight[weight_gt_0]
            _y = _y[weight_gt_0, :]
            X = X[weight_gt_0, :]
            _y = _y * np.expand_dims(sample_weight, 1)

        # we need some sort of &#34;offset&#34; array for glmnet
        # an array of shape (n_examples, n_classes)
        offset = np.zeros((X.shape[0], n_classes), dtype=np.float64, order=&#34;F&#34;)

        # You should have thought of that before you got here.
        exclude_vars = 0

        # how much each feature should be penalized relative to the others
        # this may be useful to expose to the caller if there are vars that
        # must be included in the final model or there is some prior knowledge
        # about how important some vars are relative to others, see the glmnet
        # vignette:
        # http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
        if relative_penalties is None:
            relative_penalties = np.ones(X.shape[1], dtype=np.float64, order=&#34;F&#34;)

        coef_bounds = np.empty((2, X.shape[1]), dtype=np.float64, order=&#34;F&#34;)
        coef_bounds[0, :] = self.lower_limits
        coef_bounds[1, :] = self.upper_limits

        if n_classes == 2:
            # binomial, tell glmnet there is only one class
            # otherwise we will get a coef matrix with two dimensions
            # where each pair are equal in magnitude and opposite in sign
            # also since the magnitudes are constrained to sum to one, the
            # returned coefficients would be one half of the proper values
            n_classes = 1

        # This is a stopping criterion (nx)
        # R defaults to nx = num_features, and ne = num_features + 1
        if self.max_features is None:
            max_features = X.shape[1]
        else:
            max_features = self.max_features

        if issparse(X):
            _x = csc_matrix(X, dtype=np.float64, copy=True)

            (
                self.n_lambda_,
                self.intercept_path_,
                ca,
                ia,
                nin,
                _,  # dev0
                _,  # dev
                self.lambda_path_,
                _,  # nlp
                jerr,
            ) = splognet(
                self.alpha,
                _x.shape[0],
                _x.shape[1],
                n_classes,
                _x.data,
                _x.indptr + 1,  # Fortran uses 1-based indexing
                _x.indices + 1,
                _y,
                offset,
                exclude_vars,
                relative_penalties,
                coef_bounds,
                max_features,
                X.shape[1] + 1,
                min_lambda_ratio,
                self.lambda_path,
                self.tol,
                n_lambda,
                self.standardize,
                self.fit_intercept,
                self.max_iter,
                0,
            )
        else:  # not sparse
            # some notes: glmnet requires both x and y to be float64, the two
            # arrays
            # may also be overwritten during the fitting process, so they need
            # to be copied prior to calling lognet. The fortran wrapper will
            # copy any arrays passed to a wrapped function if they are not in
            # the fortran layout, to avoid making extra copies, ensure x and y
            # are `F_CONTIGUOUS` prior to calling lognet.
            _x = X.astype(dtype=np.float64, order=&#34;F&#34;, copy=True)

            (
                self.n_lambda_,
                self.intercept_path_,
                ca,
                ia,
                nin,
                _,  # dev0
                _,  # dev
                self.lambda_path_,
                _,  # nlp
                jerr,
            ) = lognet(
                self.alpha,
                n_classes,
                _x,
                _y,
                offset,
                exclude_vars,
                relative_penalties,
                coef_bounds,
                X.shape[1] + 1,
                min_lambda_ratio,
                self.lambda_path,
                self.tol,
                max_features,
                n_lambda,
                self.standardize,
                self.fit_intercept,
                self.max_iter,
                0,
            )

        # raises RuntimeError if self.jerr_ is nonzero
        self.jerr_ = jerr
        _check_error_flag(self.jerr_, verbose_convergence=self.verbose)

        # glmnet may not return the requested number of lambda values, so we
        # need to trim the trailing zeros from the returned path so
        # len(lambda_path_) is equal to n_lambda_
        self.lambda_path_ = self.lambda_path_[: self.n_lambda_]
        # also fix the first value of lambda
        self.lambda_path_ = _fix_lambda_path(self.lambda_path_)
        self.intercept_path_ = self.intercept_path_[:, : self.n_lambda_]
        # also trim the compressed coefficient matrix
        ca = ca[:, :, : self.n_lambda_]
        # and trim the array of n_coef per lambda (may or may not be non-zero)
        nin = nin[: self.n_lambda_]
        # decompress the coefficients returned by glmnet, see doc.py
        self.coef_path_ = lsolns(X.shape[1], ca, ia, nin)
        # coef_path_ has shape (n_features, n_classes, n_lambda), we should
        # match shape for scikit-learn models:
        # (n_classes, n_features, n_lambda)
        self.coef_path_ = np.transpose(self.coef_path_, axes=(1, 0, 2))

        return self

    def decision_function(self, X, lamb=None):
        lambda_best = None
        if hasattr(self, &#34;lambda_best_&#34;):
            lambda_best = self.lambda_best_

        lamb = _check_user_lambda(self.lambda_path_, lambda_best, lamb)
        coef, intercept = _interpolate_model(
            self.lambda_path_, self.coef_path_, self.intercept_path_, lamb
        )

        # coef must be (n_classes, n_features, n_lambda)
        if coef.ndim != 3:
            # we must be working with an intercept only model
            coef = coef[:, :, np.newaxis]
        # intercept must be (n_classes, n_lambda)
        if intercept.ndim != 2:
            intercept = intercept[:, np.newaxis]

        X = check_array(X, accept_sparse=&#34;csr&#34;)
        # return (n_samples, n_classes, n_lambda)
        z = np.empty((X.shape[0], coef.shape[0], coef.shape[-1]))
        # well... sometimes we just need a for loop
        for c in range(coef.shape[0]):  # all classes
            for l in range(coef.shape[-1]):  # all values of lambda
                z[:, c, l] = X.dot(coef[c, :, l])
        z += intercept

        # drop the last dimension (lambda) when we are predicting for a single
        # value of lambda, and drop the middle dimension (class) when we are
        # predicting from a binomial model (for consistency with scikit-learn)
        return z.squeeze()

    def predict_proba(self, X, lamb=None):
        &#34;&#34;&#34;Probability estimates for each class given X.

        The returned estimates are in the same order as the values in
        classes_.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)

        lamb : array, shape (n_lambda,)
            Values of lambda from lambda_path_ from which to make predictions.
            If no values are provided, the returned predictions will be those
            corresponding to lambda_best_. The values of lamb must also be in
            the range of lambda_path_, values greater than max(lambda_path_)
            or less than  min(lambda_path_) will be clipped.

        Returns
        -------
        T : array, shape (n_samples, n_classes) or (n_samples, n_classes, n_lambda)
        &#34;&#34;&#34;
        z = self.decision_function(X, lamb)
        expit(z, z)

        # reshape z to (n_samples, n_classes, n_lambda)
        n_lambda = len(np.atleast_1d(lamb))
        z = z.reshape(X.shape[0], -1, n_lambda)

        if z.shape[1] == 1:
            # binomial, for consistency and to match scikit-learn, add the
            # complement so z has shape (n_samples, 2, n_lambda)
            z = np.concatenate((1 - z, z), axis=1)
        else:
            # normalize for multinomial
            z /= np.expand_dims(z.sum(axis=1), axis=1)

        if n_lambda == 1:
            z = z.squeeze(axis=-1)
        return z

    def predict(self, X, lamb=None):
        &#34;&#34;&#34;Predict class labels for samples in X.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)

        lamb : array, shape (n_lambda,)
            Values of lambda from lambda_path_ from which to make predictions.
            If no values are provided for lamb, the returned predictions will
            be those corresponding to lambda_best_. The values of lamb must
            also be in the range of lambda_path_, values greater than
            max(lambda_path_) or less than  min(lambda_path_) will be clipped.

        Returns
        -------
        T : array, shape (n_samples,) or (n_samples, n_lambda)
            Predicted class labels for each sample given each value of lambda
        &#34;&#34;&#34;

        scores = self.predict_proba(X, lamb)
        indices = scores.argmax(axis=1)

        return self.classes_[indices]

    def score(self, X, y, lamb=None):
        &#34;&#34;&#34;Returns the mean accuracy on the given test data and labels.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            Test samples

        y : array, shape (n_samples,)
            True labels for X

        lamb : array, shape (n_lambda,)
            Values from lambda_path_ for which to score predictions.

        Returns
        -------
        score : array, shape (n_lambda,)
            Mean accuracy for each value of lambda.
        &#34;&#34;&#34;
        pred = self.predict(X, lamb=lamb)
        return np.apply_along_axis(accuracy_score, 0, pred, y)</code></pre>
</details>
</section>
<section>
</section>
<section>
</section>
<section>
</section>
<section>
<h2 class="section-title" id="header-classes">Classes</h2>
<dl>
<dt id="dropkick.logistic.LogitNet"><code class="flex name class">
<span>class <span class="ident">LogitNet</span></span>
<span>(</span><span>alpha=1, n_lambda=100, min_lambda_ratio=0.0001, lambda_path=None, standardize=True, fit_intercept=True, lower_limits=-inf, upper_limits=inf, cut_point=1.0, n_splits=3, scoring=None, n_jobs=1, tol=1e-07, max_iter=100000, random_state=None, max_features=None, verbose=False)</span>
</code></dt>
<dd>
<div class="desc"><p>Logistic Regression with elastic net penalty.
This is a wrapper for the glmnet function lognet.</p>
<h2 id="parameters">Parameters</h2>
<dl>
<dt><strong><code>alpha</code></strong> :&ensp;<code>float</code>, default <code>1</code></dt>
<dd>The alpha parameter, 0 &lt;= alpha &lt;= 1, 0 for ridge, 1 for lasso</dd>
<dt><strong><code>n_lambda</code></strong> :&ensp;<code>int</code>, default <code>100</code></dt>
<dd>Maximum number of lambda values to compute</dd>
<dt><strong><code>min_lambda_ratio</code></strong> :&ensp;<code>float</code>, default <code>1e-4</code></dt>
<dd>In combination with <code>n_lambda</code>, the ratio of the smallest and largest
values of lambda computed.</dd>
<dt><strong><code>lambda_path</code></strong> :&ensp;<code>array</code>, default <code>None</code></dt>
<dd>In place of supplying n_lambda, provide an array of specific values
to compute. The specified values must be in decreasing order. When
None, the path of lambda values will be determined automatically. A
maximum of <code>n_lambda</code> values will be computed.</dd>
<dt><strong><code>standardize</code></strong> :&ensp;<code>bool</code>, default <code>True</code></dt>
<dd>Standardize input features prior to fitting. The final coefficients
will be on the scale of the original data regardless of the value
of standardize.</dd>
<dt><strong><code>fit_intercept</code></strong> :&ensp;<code>bool</code>, default <code>True</code></dt>
<dd>Include an intercept term in the model.</dd>
<dt><strong><code>lower_limits</code></strong> :&ensp;<code>array, (shape n_features,) default -infinity</code></dt>
<dd>Array of lower limits for each coefficient, must be non-positive.
Can be a single value (which is then replicated), else an array
corresponding to the number of features.</dd>
<dt><strong><code>upper_limits</code></strong> :&ensp;<code>array, (shape n_features,) default +infinity</code></dt>
<dd>Array of upper limits for each coefficient, must be positive.
See lower_limits.</dd>
<dt><strong><code>cut_point</code></strong> :&ensp;<code>float</code>, default <code>1</code></dt>
<dd>The cut point to use for selecting lambda_best.
arg_max lambda
cv_score(lambda) &gt;= cv_score(lambda_max) - cut_point * standard_error(lambda_max)</dd>
<dt><strong><code>n_splits</code></strong> :&ensp;<code>int</code>, default <code>3</code></dt>
<dd>Number of cross validation folds for computing performance metrics and
determining <code>lambda_best_</code> and <code>lambda_max_</code>. If non-zero, must be
at least 3.</dd>
<dt><strong><code>scoring</code></strong> :&ensp;<code>string, callable</code> or <code>None</code>, default <code>None</code></dt>
<dd>Scoring method for model selection during cross validation. When None,
defaults to classification score. Valid options are <code>accuracy</code>,
<code>roc_auc</code>, <code>average_precision</code>, <code>precision</code>, <code>recall</code>. Alternatively,
supply a function or callable object with the following signature
<code>scorer(estimator, X, y)</code>. Note, the scoring function affects the
selection of <code>lambda_best_</code> and <code>lambda_max_</code>, fitting the same data
with different scoring methods will result in the selection of
different models.</dd>
<dt><strong><code>n_jobs</code></strong> :&ensp;<code>int</code>, default <code>1</code></dt>
<dd>Maximum number of threads for computing cross validation metrics.</dd>
<dt><strong><code>tol</code></strong> :&ensp;<code>float</code>, default <code>1e-7</code></dt>
<dd>Convergence tolerance.</dd>
<dt><strong><code>max_iter</code></strong> :&ensp;<code>int</code>, default <code>100000</code></dt>
<dd>Maximum passes over the data.</dd>
<dt><strong><code>random_state</code></strong> :&ensp;<code>number</code>, default <code>None</code></dt>
<dd>Seed for the random number generator. The glmnet solver is not
deterministic, this seed is used for determining the cv folds.</dd>
<dt><strong><code>max_features</code></strong> :&ensp;<code>int</code></dt>
<dd>Optional maximum number of features with nonzero coefficients after
regularization. If not set, defaults to X.shape[1] during fit
Note, this will be ignored if the user specifies lambda_path.</dd>
<dt><strong><code>verbose</code></strong> :&ensp;<code>bool</code>, default <code>False</code></dt>
<dd>When True some warnings and log messages are suppressed.</dd>
</dl>
<h2 id="attributes">Attributes</h2>
<dl>
<dt><strong><code>classes_</code></strong> :&ensp;<code>array, shape(n_classes,)</code></dt>
<dd>The distinct classes/labels found in y.</dd>
<dt><strong><code>n_lambda_</code></strong> :&ensp;<code>int</code></dt>
<dd>The number of lambda values found by glmnet. Note, this may be less
than the number specified via n_lambda.</dd>
<dt><strong><code>lambda_path_</code></strong> :&ensp;<code>array, shape (n_lambda_,)</code></dt>
<dd>The values of lambda found by glmnet, in decreasing order.</dd>
<dt><strong><code>coef_path_</code></strong> :&ensp;<code>array, shape (n_classes, n_features, n_lambda_)</code></dt>
<dd>The set of coefficients for each value of lambda in lambda_path_.</dd>
<dt><strong><code>coef_</code></strong> :&ensp;<code>array, shape (n_clases, n_features)</code></dt>
<dd>The coefficients corresponding to lambda_best_.</dd>
<dt><strong><code>intercept_</code></strong> :&ensp;<code>array, shape (n_classes,)</code></dt>
<dd>The intercept corresponding to lambda_best_.</dd>
<dt><strong><code>intercept_path_</code></strong> :&ensp;<code>array, shape (n_classes, n_lambda_)</code></dt>
<dd>The set of intercepts for each value of lambda in lambda_path_.</dd>
<dt><strong><code>cv_mean_score_</code></strong> :&ensp;<code>array, shape (n_lambda_,)</code></dt>
<dd>The mean cv score for each value of lambda. This will be set by fit_cv.</dd>
<dt><strong><code>cv_standard_error_</code></strong> :&ensp;<code>array, shape (n_lambda_,)</code></dt>
<dd>The standard error of the mean cv score for each value of lambda, this
will be set by fit_cv.</dd>
<dt><strong><code>lambda_max_</code></strong> :&ensp;<code>float</code></dt>
<dd>The value of lambda that gives the best performance in cross
validation.</dd>
<dt><strong><code>lambda_best_</code></strong> :&ensp;<code>float</code></dt>
<dd>The largest value of lambda which is greater than lambda_max_ and
performs within cut_point * standard error of lambda_max_.</dd>
</dl></div>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">class LogitNet(BaseEstimator):
    &#34;&#34;&#34;
    Logistic Regression with elastic net penalty.
    This is a wrapper for the glmnet function lognet.

    Parameters
    ----------
    alpha : float, default 1
        The alpha parameter, 0 &lt;= alpha &lt;= 1, 0 for ridge, 1 for lasso

    n_lambda : int, default 100
        Maximum number of lambda values to compute

    min_lambda_ratio : float, default 1e-4
        In combination with `n_lambda`, the ratio of the smallest and largest
        values of lambda computed.

    lambda_path : array, default None
        In place of supplying n_lambda, provide an array of specific values
        to compute. The specified values must be in decreasing order. When
        None, the path of lambda values will be determined automatically. A
        maximum of `n_lambda` values will be computed.

    standardize : bool, default True
        Standardize input features prior to fitting. The final coefficients
        will be on the scale of the original data regardless of the value
        of standardize.

    fit_intercept : bool, default True
        Include an intercept term in the model.

    lower_limits : array, (shape n_features,) default -infinity
        Array of lower limits for each coefficient, must be non-positive.
        Can be a single value (which is then replicated), else an array
        corresponding to the number of features.

    upper_limits : array, (shape n_features,) default +infinity
        Array of upper limits for each coefficient, must be positive.
        See lower_limits.

    cut_point : float, default 1
        The cut point to use for selecting lambda_best.
            arg_max lambda  cv_score(lambda) &gt;= cv_score(lambda_max) - cut_point * standard_error(lambda_max)

    n_splits : int, default 3
        Number of cross validation folds for computing performance metrics and
        determining `lambda_best_` and `lambda_max_`. If non-zero, must be
        at least 3.

    scoring : string, callable or None, default None
        Scoring method for model selection during cross validation. When None,
        defaults to classification score. Valid options are `accuracy`,
        `roc_auc`, `average_precision`, `precision`, `recall`. Alternatively,
        supply a function or callable object with the following signature
        ``scorer(estimator, X, y)``. Note, the scoring function affects the
        selection of `lambda_best_` and `lambda_max_`, fitting the same data
        with different scoring methods will result in the selection of
        different models.

    n_jobs : int, default 1
        Maximum number of threads for computing cross validation metrics.

    tol : float, default 1e-7
        Convergence tolerance.

    max_iter : int, default 100000
        Maximum passes over the data.

    random_state : number, default None
        Seed for the random number generator. The glmnet solver is not
        deterministic, this seed is used for determining the cv folds.

    max_features : int
        Optional maximum number of features with nonzero coefficients after
        regularization. If not set, defaults to X.shape[1] during fit
        Note, this will be ignored if the user specifies lambda_path.

    verbose : bool, default False
        When True some warnings and log messages are suppressed.

    Attributes
    ----------
    classes_ : array, shape(n_classes,)
        The distinct classes/labels found in y.

    n_lambda_ : int
        The number of lambda values found by glmnet. Note, this may be less
        than the number specified via n_lambda.

    lambda_path_ : array, shape (n_lambda_,)
        The values of lambda found by glmnet, in decreasing order.

    coef_path_ : array, shape (n_classes, n_features, n_lambda_)
        The set of coefficients for each value of lambda in lambda_path_.

    coef_ : array, shape (n_clases, n_features)
        The coefficients corresponding to lambda_best_.

    intercept_ : array, shape (n_classes,)
        The intercept corresponding to lambda_best_.

    intercept_path_ : array, shape (n_classes, n_lambda_)
        The set of intercepts for each value of lambda in lambda_path_.

    cv_mean_score_ : array, shape (n_lambda_,)
        The mean cv score for each value of lambda. This will be set by fit_cv.

    cv_standard_error_ : array, shape (n_lambda_,)
        The standard error of the mean cv score for each value of lambda, this
        will be set by fit_cv.

    lambda_max_ : float
        The value of lambda that gives the best performance in cross
        validation.

    lambda_best_ : float
        The largest value of lambda which is greater than lambda_max_ and
        performs within cut_point * standard error of lambda_max_.
    &#34;&#34;&#34;

    CV = StratifiedKFold

    def __init__(
        self,
        alpha=1,
        n_lambda=100,
        min_lambda_ratio=1e-4,
        lambda_path=None,
        standardize=True,
        fit_intercept=True,
        lower_limits=-np.inf,
        upper_limits=np.inf,
        cut_point=1.0,
        n_splits=3,
        scoring=None,
        n_jobs=1,
        tol=1e-7,
        max_iter=100000,
        random_state=None,
        max_features=None,
        verbose=False,
    ):

        self.alpha = alpha
        self.n_lambda = n_lambda
        self.min_lambda_ratio = min_lambda_ratio
        self.lambda_path = lambda_path
        self.standardize = standardize
        self.lower_limits = lower_limits
        self.upper_limits = upper_limits
        self.fit_intercept = fit_intercept
        self.cut_point = cut_point
        self.n_splits = n_splits
        self.scoring = scoring
        self.n_jobs = n_jobs
        self.tol = tol
        self.max_iter = max_iter
        self.random_state = random_state
        self.max_features = max_features
        self.verbose = verbose

    def fit(self, adata, y, n_hvgs, sample_weight=None, relative_penalties=None):
        &#34;&#34;&#34;
        Fit the model to training data. If n_splits &gt; 1 also run n-fold cross
        validation on all values in lambda_path.

        The model will be fit n+1 times. On the first pass, the lambda_path
        will be determined, on the remaining passes, the model performance for
        each value of lambda. After cross validation, the attribute
        `cv_mean_score_` will contain the mean score over all folds for each
        value of lambda, and `cv_standard_error_` will contain the standard
        error of `cv_mean_score_` for each value of lambda. The value of lambda
        which achieves the best performance in cross validation will be saved
        to `lambda_max_` additionally, the largest value of lambda s.t.:
            cv_score(l) &gt;= cv_score(lambda_max_) -\
                           cut_point * standard_error(lambda_max_)
        will be saved to `lambda_best_`.

        Parameters
        ----------
        adata : anndata.AnnData, shape (n_samples, n_features)
            Input features in .X

        y : array, shape (n_samples,)
            Target values

        n_hvgs : int, number of highly-variable genes to feature-select

        sample_weight : array, shape (n_samples,)
            Optional weight vector for observations

        relative_penalties: array, shape (n_features,)
            Optional relative weight vector for penalty.
            0 entries remove penalty.

        Returns
        -------
        self : object
            Returns self.
        &#34;&#34;&#34;
        # X is scaled counts for all cells in HVGs only for first run
        X = adata[:, adata.var.highly_variable].X.copy()

        X, y = check_X_y(X, y, accept_sparse=&#34;csr&#34;, ensure_min_samples=2)
        if sample_weight is None:
            sample_weight = np.ones(X.shape[0])
        else:
            sample_weight = np.asarray(sample_weight)

        if not np.isscalar(self.lower_limits):
            self.lower_limits = np.asarray(self.lower_limits)
            if len(self.lower_limits) != X.shape[1]:
                raise ValueError(&#34;lower_limits must equal number of features&#34;)

        if not np.isscalar(self.upper_limits):
            self.upper_limits = np.asarray(self.upper_limits)
            if len(self.upper_limits) != X.shape[1]:
                raise ValueError(&#34;upper_limits must equal number of features&#34;)

        if (
            any(self.lower_limits &gt; 0)
            if isinstance(self.lower_limits, np.ndarray)
            else self.lower_limits &gt; 0
        ):
            raise ValueError(&#34;lower_limits must be non-positive&#34;)

        if (
            any(self.upper_limits &lt; 0)
            if isinstance(self.upper_limits, np.ndarray)
            else self.upper_limits &lt; 0
        ):
            raise ValueError(&#34;upper_limits must be positive&#34;)

        if self.alpha &gt; 1 or self.alpha &lt; 0:
            raise ValueError(&#34;alpha must be between 0 and 1&#34;)

        # fit the model
        self._fit(X, y, sample_weight, relative_penalties)

        # score each model on the path of lambda values found by glmnet and
        # select the best scoring
        if self.n_splits &gt;= 3:
            self._cv = self.CV(
                n_splits=self.n_splits, shuffle=True, random_state=self.random_state
            )

            cv_scores, _hvgs = _score_lambda_path(
                self,
                adata,
                y,
                n_hvgs,
                sample_weight,
                relative_penalties,
                self.scoring,
                n_jobs=self.n_jobs,
                verbose=self.verbose,
            )

            self.cv_mean_score_ = np.atleast_1d(np.mean(cv_scores, axis=0))
            self.cv_standard_error_ = np.atleast_1d(stats.sem(cv_scores))

            self.lambda_max_inx_ = np.argmax(self.cv_mean_score_)
            self.lambda_max_ = self.lambda_path_[self.lambda_max_inx_]

            target_score = (
                self.cv_mean_score_[self.lambda_max_inx_]
                - self.cut_point * self.cv_standard_error_[self.lambda_max_inx_]
            )

            self.lambda_best_inx_ = np.argwhere(self.cv_mean_score_ &gt;= target_score)[0]
            self.lambda_best_ = self.lambda_path_[self.lambda_best_inx_]

            self.coef_ = self.coef_path_[..., self.lambda_best_inx_]
            self.coef_ = self.coef_.squeeze(axis=self.coef_.ndim - 1)
            self.intercept_ = self.intercept_path_[..., self.lambda_best_inx_].squeeze()
            if self.intercept_.shape == ():  # convert 0d array to scalar
                self.intercept_ = float(self.intercept_)

        return self

    def _fit(self, X, y, sample_weight=None, relative_penalties=None):
        if self.lambda_path is not None:
            n_lambda = len(self.lambda_path)
            min_lambda_ratio = 1.0
        else:
            n_lambda = self.n_lambda
            min_lambda_ratio = self.min_lambda_ratio

        check_classification_targets(y)
        self.classes_ = np.unique(y)  # the output of np.unique is sorted
        n_classes = len(self.classes_)
        if n_classes &lt; 2:
            raise ValueError(&#34;Training data need to contain at least 2 &#34; &#34;classes.&#34;)

        # glmnet requires the labels a one-hot-encoded array of
        # (n_samples, n_classes)
        if n_classes == 2:
            # Normally we use 1/0 for the positive and negative classes. Since
            # np.unique sorts the output, the negative class will be in the 0th
            # column. We want a model predicting the positive class, not the
            # negative class, so we flip the columns here (the != condition).
            #
            # Broadcast comparison of self.classes_ to all rows of y. See the
            # numpy rules on broadcasting for more info, essentially this
            # &#34;reshapes&#34; y to (n_samples, n_classes) and self.classes_ to
            # (n_samples, n_classes) and performs an element-wise comparison
            # resulting in _y with shape (n_samples, n_classes).
            _y = (y[:, None] != self.classes_).astype(np.float64, order=&#34;F&#34;)
        else:
            # multinomial case, glmnet uses the entire array so we can
            # keep the original order.
            _y = (y[:, None] == self.classes_).astype(np.float64, order=&#34;F&#34;)

        # use sample weights, making sure all weights are positive
        # this is inspired by the R wrapper for glmnet, in lognet.R
        if sample_weight is not None:
            weight_gt_0 = sample_weight &gt; 0
            sample_weight = sample_weight[weight_gt_0]
            _y = _y[weight_gt_0, :]
            X = X[weight_gt_0, :]
            _y = _y * np.expand_dims(sample_weight, 1)

        # we need some sort of &#34;offset&#34; array for glmnet
        # an array of shape (n_examples, n_classes)
        offset = np.zeros((X.shape[0], n_classes), dtype=np.float64, order=&#34;F&#34;)

        # You should have thought of that before you got here.
        exclude_vars = 0

        # how much each feature should be penalized relative to the others
        # this may be useful to expose to the caller if there are vars that
        # must be included in the final model or there is some prior knowledge
        # about how important some vars are relative to others, see the glmnet
        # vignette:
        # http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
        if relative_penalties is None:
            relative_penalties = np.ones(X.shape[1], dtype=np.float64, order=&#34;F&#34;)

        coef_bounds = np.empty((2, X.shape[1]), dtype=np.float64, order=&#34;F&#34;)
        coef_bounds[0, :] = self.lower_limits
        coef_bounds[1, :] = self.upper_limits

        if n_classes == 2:
            # binomial, tell glmnet there is only one class
            # otherwise we will get a coef matrix with two dimensions
            # where each pair are equal in magnitude and opposite in sign
            # also since the magnitudes are constrained to sum to one, the
            # returned coefficients would be one half of the proper values
            n_classes = 1

        # This is a stopping criterion (nx)
        # R defaults to nx = num_features, and ne = num_features + 1
        if self.max_features is None:
            max_features = X.shape[1]
        else:
            max_features = self.max_features

        if issparse(X):
            _x = csc_matrix(X, dtype=np.float64, copy=True)

            (
                self.n_lambda_,
                self.intercept_path_,
                ca,
                ia,
                nin,
                _,  # dev0
                _,  # dev
                self.lambda_path_,
                _,  # nlp
                jerr,
            ) = splognet(
                self.alpha,
                _x.shape[0],
                _x.shape[1],
                n_classes,
                _x.data,
                _x.indptr + 1,  # Fortran uses 1-based indexing
                _x.indices + 1,
                _y,
                offset,
                exclude_vars,
                relative_penalties,
                coef_bounds,
                max_features,
                X.shape[1] + 1,
                min_lambda_ratio,
                self.lambda_path,
                self.tol,
                n_lambda,
                self.standardize,
                self.fit_intercept,
                self.max_iter,
                0,
            )
        else:  # not sparse
            # some notes: glmnet requires both x and y to be float64, the two
            # arrays
            # may also be overwritten during the fitting process, so they need
            # to be copied prior to calling lognet. The fortran wrapper will
            # copy any arrays passed to a wrapped function if they are not in
            # the fortran layout, to avoid making extra copies, ensure x and y
            # are `F_CONTIGUOUS` prior to calling lognet.
            _x = X.astype(dtype=np.float64, order=&#34;F&#34;, copy=True)

            (
                self.n_lambda_,
                self.intercept_path_,
                ca,
                ia,
                nin,
                _,  # dev0
                _,  # dev
                self.lambda_path_,
                _,  # nlp
                jerr,
            ) = lognet(
                self.alpha,
                n_classes,
                _x,
                _y,
                offset,
                exclude_vars,
                relative_penalties,
                coef_bounds,
                X.shape[1] + 1,
                min_lambda_ratio,
                self.lambda_path,
                self.tol,
                max_features,
                n_lambda,
                self.standardize,
                self.fit_intercept,
                self.max_iter,
                0,
            )

        # raises RuntimeError if self.jerr_ is nonzero
        self.jerr_ = jerr
        _check_error_flag(self.jerr_, verbose_convergence=self.verbose)

        # glmnet may not return the requested number of lambda values, so we
        # need to trim the trailing zeros from the returned path so
        # len(lambda_path_) is equal to n_lambda_
        self.lambda_path_ = self.lambda_path_[: self.n_lambda_]
        # also fix the first value of lambda
        self.lambda_path_ = _fix_lambda_path(self.lambda_path_)
        self.intercept_path_ = self.intercept_path_[:, : self.n_lambda_]
        # also trim the compressed coefficient matrix
        ca = ca[:, :, : self.n_lambda_]
        # and trim the array of n_coef per lambda (may or may not be non-zero)
        nin = nin[: self.n_lambda_]
        # decompress the coefficients returned by glmnet, see doc.py
        self.coef_path_ = lsolns(X.shape[1], ca, ia, nin)
        # coef_path_ has shape (n_features, n_classes, n_lambda), we should
        # match shape for scikit-learn models:
        # (n_classes, n_features, n_lambda)
        self.coef_path_ = np.transpose(self.coef_path_, axes=(1, 0, 2))

        return self

    def decision_function(self, X, lamb=None):
        lambda_best = None
        if hasattr(self, &#34;lambda_best_&#34;):
            lambda_best = self.lambda_best_

        lamb = _check_user_lambda(self.lambda_path_, lambda_best, lamb)
        coef, intercept = _interpolate_model(
            self.lambda_path_, self.coef_path_, self.intercept_path_, lamb
        )

        # coef must be (n_classes, n_features, n_lambda)
        if coef.ndim != 3:
            # we must be working with an intercept only model
            coef = coef[:, :, np.newaxis]
        # intercept must be (n_classes, n_lambda)
        if intercept.ndim != 2:
            intercept = intercept[:, np.newaxis]

        X = check_array(X, accept_sparse=&#34;csr&#34;)
        # return (n_samples, n_classes, n_lambda)
        z = np.empty((X.shape[0], coef.shape[0], coef.shape[-1]))
        # well... sometimes we just need a for loop
        for c in range(coef.shape[0]):  # all classes
            for l in range(coef.shape[-1]):  # all values of lambda
                z[:, c, l] = X.dot(coef[c, :, l])
        z += intercept

        # drop the last dimension (lambda) when we are predicting for a single
        # value of lambda, and drop the middle dimension (class) when we are
        # predicting from a binomial model (for consistency with scikit-learn)
        return z.squeeze()

    def predict_proba(self, X, lamb=None):
        &#34;&#34;&#34;Probability estimates for each class given X.

        The returned estimates are in the same order as the values in
        classes_.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)

        lamb : array, shape (n_lambda,)
            Values of lambda from lambda_path_ from which to make predictions.
            If no values are provided, the returned predictions will be those
            corresponding to lambda_best_. The values of lamb must also be in
            the range of lambda_path_, values greater than max(lambda_path_)
            or less than  min(lambda_path_) will be clipped.

        Returns
        -------
        T : array, shape (n_samples, n_classes) or (n_samples, n_classes, n_lambda)
        &#34;&#34;&#34;
        z = self.decision_function(X, lamb)
        expit(z, z)

        # reshape z to (n_samples, n_classes, n_lambda)
        n_lambda = len(np.atleast_1d(lamb))
        z = z.reshape(X.shape[0], -1, n_lambda)

        if z.shape[1] == 1:
            # binomial, for consistency and to match scikit-learn, add the
            # complement so z has shape (n_samples, 2, n_lambda)
            z = np.concatenate((1 - z, z), axis=1)
        else:
            # normalize for multinomial
            z /= np.expand_dims(z.sum(axis=1), axis=1)

        if n_lambda == 1:
            z = z.squeeze(axis=-1)
        return z

    def predict(self, X, lamb=None):
        &#34;&#34;&#34;Predict class labels for samples in X.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)

        lamb : array, shape (n_lambda,)
            Values of lambda from lambda_path_ from which to make predictions.
            If no values are provided for lamb, the returned predictions will
            be those corresponding to lambda_best_. The values of lamb must
            also be in the range of lambda_path_, values greater than
            max(lambda_path_) or less than  min(lambda_path_) will be clipped.

        Returns
        -------
        T : array, shape (n_samples,) or (n_samples, n_lambda)
            Predicted class labels for each sample given each value of lambda
        &#34;&#34;&#34;

        scores = self.predict_proba(X, lamb)
        indices = scores.argmax(axis=1)

        return self.classes_[indices]

    def score(self, X, y, lamb=None):
        &#34;&#34;&#34;Returns the mean accuracy on the given test data and labels.

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            Test samples

        y : array, shape (n_samples,)
            True labels for X

        lamb : array, shape (n_lambda,)
            Values from lambda_path_ for which to score predictions.

        Returns
        -------
        score : array, shape (n_lambda,)
            Mean accuracy for each value of lambda.
        &#34;&#34;&#34;
        pred = self.predict(X, lamb=lamb)
        return np.apply_along_axis(accuracy_score, 0, pred, y)</code></pre>
</details>
<h3>Ancestors</h3>
<ul class="hlist">
<li>sklearn.base.BaseEstimator</li>
</ul>
<h3>Class variables</h3>
<dl>
<dt id="dropkick.logistic.LogitNet.CV"><code class="name">var <span class="ident">CV</span></code></dt>
<dd>
<div class="desc"><p>Stratified K-Folds cross-validator.</p>
<p>Provides train/test indices to split data in train/test sets.</p>
<p>This cross-validation object is a variation of KFold that returns
stratified folds. The folds are made by preserving the percentage of
samples for each class.</p>
<p>Read more in the :ref:<code>User Guide &lt;stratified_k_fold&gt;</code>.</p>
<h2 id="parameters">Parameters</h2>
<dl>
<dt><strong><code>n_splits</code></strong> :&ensp;<code>int</code>, default=<code>5</code></dt>
<dd>
<p>Number of folds. Must be at least 2.</p>
<div class="admonition versionchanged">
<p class="admonition-title">Changed in version:&ensp;0.22</p>
<p><code>n_splits</code> default value changed from 3 to 5.</p>
</div>
</dd>
<dt><strong><code>shuffle</code></strong> :&ensp;<code>bool</code>, default=<code>False</code></dt>
<dd>Whether to shuffle each class's samples before splitting into batches.
Note that the samples within each split will not be shuffled.</dd>
<dt><strong><code>random_state</code></strong> :&ensp;<code>int, RandomState instance</code> or <code>None</code>, default=<code>None</code></dt>
<dd>When <code>shuffle</code> is True, <code>random_state</code> affects the ordering of the
indices, which controls the randomness of each fold for each class.
Otherwise, leave <code>random_state</code> as <code>None</code>.
Pass an int for reproducible output across multiple function calls.
See :term:<code>Glossary &lt;random_state&gt;</code>.</dd>
</dl>
<h2 id="examples">Examples</h2>
<pre><code class="language-python-repl">&gt;&gt;&gt; import numpy as np
&gt;&gt;&gt; from sklearn.model_selection import StratifiedKFold
&gt;&gt;&gt; X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
&gt;&gt;&gt; y = np.array([0, 0, 1, 1])
&gt;&gt;&gt; skf = StratifiedKFold(n_splits=2)
&gt;&gt;&gt; skf.get_n_splits(X, y)
2
&gt;&gt;&gt; print(skf)
StratifiedKFold(n_splits=2, random_state=None, shuffle=False)
&gt;&gt;&gt; for train_index, test_index in skf.split(X, y):
...     print(&quot;TRAIN:&quot;, train_index, &quot;TEST:&quot;, test_index)
...     X_train, X_test = X[train_index], X[test_index]
...     y_train, y_test = y[train_index], y[test_index]
TRAIN: [1 3] TEST: [0 2]
TRAIN: [0 2] TEST: [1 3]
</code></pre>
<h2 id="notes">Notes</h2>
<p>The implementation is designed to:</p>
<ul>
<li>Generate test sets such that all contain the same distribution of
classes, or as close as possible.</li>
<li>Be invariant to class label: relabelling <code>y = ["Happy", "Sad"]</code> to
<code>y = [1, 0]</code> should not change the indices generated.</li>
<li>Preserve order dependencies in the dataset ordering, when
<code>shuffle=False</code>: all samples from class k in some test set were
contiguous in y, or separated in y by samples from classes other than k.</li>
<li>Generate test sets where the smallest and largest differ by at most one
sample.</li>
</ul>
<div class="admonition versionchanged">
<p class="admonition-title">Changed in version:&ensp;0.22</p>
<p>The previous implementation did not follow the last constraint.</p>
</div>
<h2 id="see-also">See Also</h2>
<dl>
<dt><code>RepeatedStratifiedKFold</code></dt>
<dd>Repeats Stratified K-Fold n times.</dd>
</dl></div>
</dd>
</dl>
<h3>Methods</h3>
<dl>
<dt id="dropkick.logistic.LogitNet.decision_function"><code class="name flex">
<span>def <span class="ident">decision_function</span></span>(<span>self, X, lamb=None)</span>
</code></dt>
<dd>
<div class="desc"></div>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">def decision_function(self, X, lamb=None):
    lambda_best = None
    if hasattr(self, &#34;lambda_best_&#34;):
        lambda_best = self.lambda_best_

    lamb = _check_user_lambda(self.lambda_path_, lambda_best, lamb)
    coef, intercept = _interpolate_model(
        self.lambda_path_, self.coef_path_, self.intercept_path_, lamb
    )

    # coef must be (n_classes, n_features, n_lambda)
    if coef.ndim != 3:
        # we must be working with an intercept only model
        coef = coef[:, :, np.newaxis]
    # intercept must be (n_classes, n_lambda)
    if intercept.ndim != 2:
        intercept = intercept[:, np.newaxis]

    X = check_array(X, accept_sparse=&#34;csr&#34;)
    # return (n_samples, n_classes, n_lambda)
    z = np.empty((X.shape[0], coef.shape[0], coef.shape[-1]))
    # well... sometimes we just need a for loop
    for c in range(coef.shape[0]):  # all classes
        for l in range(coef.shape[-1]):  # all values of lambda
            z[:, c, l] = X.dot(coef[c, :, l])
    z += intercept

    # drop the last dimension (lambda) when we are predicting for a single
    # value of lambda, and drop the middle dimension (class) when we are
    # predicting from a binomial model (for consistency with scikit-learn)
    return z.squeeze()</code></pre>
</details>
</dd>
<dt id="dropkick.logistic.LogitNet.fit"><code class="name flex">
<span>def <span class="ident">fit</span></span>(<span>self, adata, y, n_hvgs, sample_weight=None, relative_penalties=None)</span>
</code></dt>
<dd>
<div class="desc"><p>Fit the model to training data. If n_splits &gt; 1 also run n-fold cross
validation on all values in lambda_path.</p>
<p>The model will be fit n+1 times. On the first pass, the lambda_path
will be determined, on the remaining passes, the model performance for
each value of lambda. After cross validation, the attribute
<code>cv_mean_score_</code> will contain the mean score over all folds for each
value of lambda, and <code>cv_standard_error_</code> will contain the standard
error of <code>cv_mean_score_</code> for each value of lambda. The value of lambda
which achieves the best performance in cross validation will be saved
to <code>lambda_max_</code> additionally, the largest value of lambda s.t.:
cv_score(l) &gt;= cv_score(lambda_max_) -
cut_point * standard_error(lambda_max_)
will be saved to <code>lambda_best_</code>.</p>
<h2 id="parameters">Parameters</h2>
<dl>
<dt><strong><code>adata</code></strong> :&ensp;<code>anndata.AnnData, shape (n_samples, n_features)</code></dt>
<dd>Input features in .X</dd>
<dt><strong><code>y</code></strong> :&ensp;<code>array, shape (n_samples,)</code></dt>
<dd>Target values</dd>
<dt><strong><code>n_hvgs</code></strong> :&ensp;<code>int, number</code> of <code>highly-variable genes to feature-select</code></dt>
<dd>&nbsp;</dd>
<dt><strong><code>sample_weight</code></strong> :&ensp;<code>array, shape (n_samples,)</code></dt>
<dd>Optional weight vector for observations</dd>
<dt><strong><code>relative_penalties</code></strong> :&ensp;<code>array, shape (n_features,)</code></dt>
<dd>Optional relative weight vector for penalty.
0 entries remove penalty.</dd>
</dl>
<h2 id="returns">Returns</h2>
<dl>
<dt><strong><code>self</code></strong> :&ensp;<code>object</code></dt>
<dd>Returns self.</dd>
</dl></div>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">def fit(self, adata, y, n_hvgs, sample_weight=None, relative_penalties=None):
    &#34;&#34;&#34;
    Fit the model to training data. If n_splits &gt; 1 also run n-fold cross
    validation on all values in lambda_path.

    The model will be fit n+1 times. On the first pass, the lambda_path
    will be determined, on the remaining passes, the model performance for
    each value of lambda. After cross validation, the attribute
    `cv_mean_score_` will contain the mean score over all folds for each
    value of lambda, and `cv_standard_error_` will contain the standard
    error of `cv_mean_score_` for each value of lambda. The value of lambda
    which achieves the best performance in cross validation will be saved
    to `lambda_max_` additionally, the largest value of lambda s.t.:
        cv_score(l) &gt;= cv_score(lambda_max_) -\
                       cut_point * standard_error(lambda_max_)
    will be saved to `lambda_best_`.

    Parameters
    ----------
    adata : anndata.AnnData, shape (n_samples, n_features)
        Input features in .X

    y : array, shape (n_samples,)
        Target values

    n_hvgs : int, number of highly-variable genes to feature-select

    sample_weight : array, shape (n_samples,)
        Optional weight vector for observations

    relative_penalties: array, shape (n_features,)
        Optional relative weight vector for penalty.
        0 entries remove penalty.

    Returns
    -------
    self : object
        Returns self.
    &#34;&#34;&#34;
    # X is scaled counts for all cells in HVGs only for first run
    X = adata[:, adata.var.highly_variable].X.copy()

    X, y = check_X_y(X, y, accept_sparse=&#34;csr&#34;, ensure_min_samples=2)
    if sample_weight is None:
        sample_weight = np.ones(X.shape[0])
    else:
        sample_weight = np.asarray(sample_weight)

    if not np.isscalar(self.lower_limits):
        self.lower_limits = np.asarray(self.lower_limits)
        if len(self.lower_limits) != X.shape[1]:
            raise ValueError(&#34;lower_limits must equal number of features&#34;)

    if not np.isscalar(self.upper_limits):
        self.upper_limits = np.asarray(self.upper_limits)
        if len(self.upper_limits) != X.shape[1]:
            raise ValueError(&#34;upper_limits must equal number of features&#34;)

    if (
        any(self.lower_limits &gt; 0)
        if isinstance(self.lower_limits, np.ndarray)
        else self.lower_limits &gt; 0
    ):
        raise ValueError(&#34;lower_limits must be non-positive&#34;)

    if (
        any(self.upper_limits &lt; 0)
        if isinstance(self.upper_limits, np.ndarray)
        else self.upper_limits &lt; 0
    ):
        raise ValueError(&#34;upper_limits must be positive&#34;)

    if self.alpha &gt; 1 or self.alpha &lt; 0:
        raise ValueError(&#34;alpha must be between 0 and 1&#34;)

    # fit the model
    self._fit(X, y, sample_weight, relative_penalties)

    # score each model on the path of lambda values found by glmnet and
    # select the best scoring
    if self.n_splits &gt;= 3:
        self._cv = self.CV(
            n_splits=self.n_splits, shuffle=True, random_state=self.random_state
        )

        cv_scores, _hvgs = _score_lambda_path(
            self,
            adata,
            y,
            n_hvgs,
            sample_weight,
            relative_penalties,
            self.scoring,
            n_jobs=self.n_jobs,
            verbose=self.verbose,
        )

        self.cv_mean_score_ = np.atleast_1d(np.mean(cv_scores, axis=0))
        self.cv_standard_error_ = np.atleast_1d(stats.sem(cv_scores))

        self.lambda_max_inx_ = np.argmax(self.cv_mean_score_)
        self.lambda_max_ = self.lambda_path_[self.lambda_max_inx_]

        target_score = (
            self.cv_mean_score_[self.lambda_max_inx_]
            - self.cut_point * self.cv_standard_error_[self.lambda_max_inx_]
        )

        self.lambda_best_inx_ = np.argwhere(self.cv_mean_score_ &gt;= target_score)[0]
        self.lambda_best_ = self.lambda_path_[self.lambda_best_inx_]

        self.coef_ = self.coef_path_[..., self.lambda_best_inx_]
        self.coef_ = self.coef_.squeeze(axis=self.coef_.ndim - 1)
        self.intercept_ = self.intercept_path_[..., self.lambda_best_inx_].squeeze()
        if self.intercept_.shape == ():  # convert 0d array to scalar
            self.intercept_ = float(self.intercept_)

    return self</code></pre>
</details>
</dd>
<dt id="dropkick.logistic.LogitNet.predict"><code class="name flex">
<span>def <span class="ident">predict</span></span>(<span>self, X, lamb=None)</span>
</code></dt>
<dd>
<div class="desc"><p>Predict class labels for samples in X.</p>
<h2 id="parameters">Parameters</h2>
<dl>
<dt><strong><code>X</code></strong> :&ensp;<code>array, shape (n_samples, n_features)</code></dt>
<dd>&nbsp;</dd>
<dt><strong><code>lamb</code></strong> :&ensp;<code>array, shape (n_lambda,)</code></dt>
<dd>Values of lambda from lambda_path_ from which to make predictions.
If no values are provided for lamb, the returned predictions will
be those corresponding to lambda_best_. The values of lamb must
also be in the range of lambda_path_, values greater than
max(lambda_path_) or less than
min(lambda_path_) will be clipped.</dd>
</dl>
<h2 id="returns">Returns</h2>
<dl>
<dt><strong><code>T</code></strong> :&ensp;<code>array, shape (n_samples,)</code> or <code>(n_samples, n_lambda)</code></dt>
<dd>Predicted class labels for each sample given each value of lambda</dd>
</dl></div>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">def predict(self, X, lamb=None):
    &#34;&#34;&#34;Predict class labels for samples in X.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)

    lamb : array, shape (n_lambda,)
        Values of lambda from lambda_path_ from which to make predictions.
        If no values are provided for lamb, the returned predictions will
        be those corresponding to lambda_best_. The values of lamb must
        also be in the range of lambda_path_, values greater than
        max(lambda_path_) or less than  min(lambda_path_) will be clipped.

    Returns
    -------
    T : array, shape (n_samples,) or (n_samples, n_lambda)
        Predicted class labels for each sample given each value of lambda
    &#34;&#34;&#34;

    scores = self.predict_proba(X, lamb)
    indices = scores.argmax(axis=1)

    return self.classes_[indices]</code></pre>
</details>
</dd>
<dt id="dropkick.logistic.LogitNet.predict_proba"><code class="name flex">
<span>def <span class="ident">predict_proba</span></span>(<span>self, X, lamb=None)</span>
</code></dt>
<dd>
<div class="desc"><p>Probability estimates for each class given X.</p>
<p>The returned estimates are in the same order as the values in
classes_.</p>
<h2 id="parameters">Parameters</h2>
<dl>
<dt><strong><code>X</code></strong> :&ensp;<code>array, shape (n_samples, n_features)</code></dt>
<dd>&nbsp;</dd>
<dt><strong><code>lamb</code></strong> :&ensp;<code>array, shape (n_lambda,)</code></dt>
<dd>Values of lambda from lambda_path_ from which to make predictions.
If no values are provided, the returned predictions will be those
corresponding to lambda_best_. The values of lamb must also be in
the range of lambda_path_, values greater than max(lambda_path_)
or less than
min(lambda_path_) will be clipped.</dd>
</dl>
<h2 id="returns">Returns</h2>
<dl>
<dt><strong><code>T</code></strong> :&ensp;<code>array, shape (n_samples, n_classes)</code> or <code>(n_samples, n_classes, n_lambda)</code></dt>
<dd>&nbsp;</dd>
</dl></div>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">def predict_proba(self, X, lamb=None):
    &#34;&#34;&#34;Probability estimates for each class given X.

    The returned estimates are in the same order as the values in
    classes_.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)

    lamb : array, shape (n_lambda,)
        Values of lambda from lambda_path_ from which to make predictions.
        If no values are provided, the returned predictions will be those
        corresponding to lambda_best_. The values of lamb must also be in
        the range of lambda_path_, values greater than max(lambda_path_)
        or less than  min(lambda_path_) will be clipped.

    Returns
    -------
    T : array, shape (n_samples, n_classes) or (n_samples, n_classes, n_lambda)
    &#34;&#34;&#34;
    z = self.decision_function(X, lamb)
    expit(z, z)

    # reshape z to (n_samples, n_classes, n_lambda)
    n_lambda = len(np.atleast_1d(lamb))
    z = z.reshape(X.shape[0], -1, n_lambda)

    if z.shape[1] == 1:
        # binomial, for consistency and to match scikit-learn, add the
        # complement so z has shape (n_samples, 2, n_lambda)
        z = np.concatenate((1 - z, z), axis=1)
    else:
        # normalize for multinomial
        z /= np.expand_dims(z.sum(axis=1), axis=1)

    if n_lambda == 1:
        z = z.squeeze(axis=-1)
    return z</code></pre>
</details>
</dd>
<dt id="dropkick.logistic.LogitNet.score"><code class="name flex">
<span>def <span class="ident">score</span></span>(<span>self, X, y, lamb=None)</span>
</code></dt>
<dd>
<div class="desc"><p>Returns the mean accuracy on the given test data and labels.</p>
<h2 id="parameters">Parameters</h2>
<dl>
<dt><strong><code>X</code></strong> :&ensp;<code>array, shape (n_samples, n_features)</code></dt>
<dd>Test samples</dd>
<dt><strong><code>y</code></strong> :&ensp;<code>array, shape (n_samples,)</code></dt>
<dd>True labels for X</dd>
<dt><strong><code>lamb</code></strong> :&ensp;<code>array, shape (n_lambda,)</code></dt>
<dd>Values from lambda_path_ for which to score predictions.</dd>
</dl>
<h2 id="returns">Returns</h2>
<dl>
<dt><strong><code>score</code></strong> :&ensp;<code>array, shape (n_lambda,)</code></dt>
<dd>Mean accuracy for each value of lambda.</dd>
</dl></div>
<details class="source">
<summary>
<span>Expand source code</span>
</summary>
<pre><code class="python">def score(self, X, y, lamb=None):
    &#34;&#34;&#34;Returns the mean accuracy on the given test data and labels.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        Test samples

    y : array, shape (n_samples,)
        True labels for X

    lamb : array, shape (n_lambda,)
        Values from lambda_path_ for which to score predictions.

    Returns
    -------
    score : array, shape (n_lambda,)
        Mean accuracy for each value of lambda.
    &#34;&#34;&#34;
    pred = self.predict(X, lamb=lamb)
    return np.apply_along_axis(accuracy_score, 0, pred, y)</code></pre>
</details>
</dd>
</dl>
</dd>
</dl>
</section>
</article>
<nav id="sidebar">
<h1>Index</h1>
<div class="toc">
<ul></ul>
</div>
<ul id="index">
<li><h3>Super-module</h3>
<ul>
<li><code><a title="dropkick" href="index.html">dropkick</a></code></li>
</ul>
</li>
<li><h3><a href="#header-classes">Classes</a></h3>
<ul>
<li>
<h4><code><a title="dropkick.logistic.LogitNet" href="#dropkick.logistic.LogitNet">LogitNet</a></code></h4>
<ul class="two-column">
<li><code><a title="dropkick.logistic.LogitNet.CV" href="#dropkick.logistic.LogitNet.CV">CV</a></code></li>
<li><code><a title="dropkick.logistic.LogitNet.decision_function" href="#dropkick.logistic.LogitNet.decision_function">decision_function</a></code></li>
<li><code><a title="dropkick.logistic.LogitNet.fit" href="#dropkick.logistic.LogitNet.fit">fit</a></code></li>
<li><code><a title="dropkick.logistic.LogitNet.predict" href="#dropkick.logistic.LogitNet.predict">predict</a></code></li>
<li><code><a title="dropkick.logistic.LogitNet.predict_proba" href="#dropkick.logistic.LogitNet.predict_proba">predict_proba</a></code></li>
<li><code><a title="dropkick.logistic.LogitNet.score" href="#dropkick.logistic.LogitNet.score">score</a></code></li>
</ul>
</li>
</ul>
</li>
</ul>
</nav>
</main>
<footer id="footer">
<p>Generated by <a href="https://pdoc3.github.io/pdoc"><cite>pdoc</cite> 0.9.2</a>.</p>
</footer>
</body>
</html>
