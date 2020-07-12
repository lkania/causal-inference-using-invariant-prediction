#### Causal inference using invariant prediction

Given a target covariate $$Y$$, and a set of covariates $$X$$. Our goal is to discover if a causal set of covariates exists in $$X$$.

**Main assumption** Let $$Y^{(i)}$$ denote the distribution of the target variable when $$X_i \in X$$ was intervened, and $$Y^{(o)}$$ denote the distribution of the target variable when no covariate was intervened. For the true causal set $$S$$ (which might not be a subset of $$X$$), the following holds

\begin{equation}
\forall_{i} \quad Y^{(o)} | S^{(o)} \equiv Y^{(i)} | S^{(i)}
\end{equation}

Furthermore, assuming an underlying Linear SEM model with additive normal noise, the assumption boils down to

\begin{equation}\label{sem}
\forall_{i} \quad \gamma^{(o)} = \gamma^{(i)}
\end{equation}

where $$\gamma^{(i)}$$ corresponds to the coefficients resulting of linearly regressing $$Y^{(i)}$$ on $$S^{(i)}$$.

An algorithm for finding $$S$$ under the above-mentioned assumptions was proposed by [Peters et al.(2016)]( https://doi.org/10.1111/rssb.12167) for the case when the experimenters already have a sample of the joint distribution of all covariates (i.e. a sample from the *observational distribution*), and samples from different interventions at one or more covariates from $$X$$. A concise reimplementation of their method and some examples with synthetic datasets are provided at [[code]](https://github.com/lkania/causal-inference-using-invariant-prediction/blob/master/causality.R)
