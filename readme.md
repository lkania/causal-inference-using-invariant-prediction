#### Causal inference using invariant prediction

Given a target covariate <img src="/tex/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode&sanitize=true" align=middle width=13.19638649999999pt height=22.465723500000017pt/>, and a set of covariates <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>. Our goal is to discover if a causal set of covariates exists in <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>.

**Main assumption** Let <img src="/tex/0218f905ddd9de68c5116265ca5d37e5.svg?invert_in_darkmode&sanitize=true" align=middle width=28.12131134999999pt height=29.190975000000005pt/> denote the distribution of the target variable when <img src="/tex/e08127ac8bf4080a89dbb978936c4bb9.svg?invert_in_darkmode&sanitize=true" align=middle width=54.091374149999986pt height=22.465723500000017pt/> was intervened, and <img src="/tex/1aa79116ffbe8e8283521cb5813e7bb1.svg?invert_in_darkmode&sanitize=true" align=middle width=29.95901534999999pt height=29.190975000000005pt/> denote the distribution of the target variable when no covariate was intervened. For the true causal set <img src="/tex/e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode&sanitize=true" align=middle width=11.027402099999989pt height=22.465723500000017pt/> (which might not be a subset of <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>), the following holds

<p align="center"><img src="/tex/a6bb26e08efe3508cf786bc6e1670073.svg?invert_in_darkmode&sanitize=true" align=middle width=438.73890059999997pt height=19.526994300000002pt/></p>

Furthermore, assuming an underlying Linear SEM model with additive normal noise, the assumption boils down to

<p align="center"><img src="/tex/5f9d537b028ccd46fd9734ac200f5ff1.svg?invert_in_darkmode&sanitize=true" align=middle width=402.70715594999996pt height=19.526994300000002pt/></p>

where <img src="/tex/97bbfccbf4c5b284f8a5948bd01bbe31.svg?invert_in_darkmode&sanitize=true" align=middle width=24.34879259999999pt height=29.190975000000005pt/> corresponds to the coefficients resulting of linearly regressing <img src="/tex/0218f905ddd9de68c5116265ca5d37e5.svg?invert_in_darkmode&sanitize=true" align=middle width=28.12131134999999pt height=29.190975000000005pt/> on <img src="/tex/886064aa073cd91c17d30ddfd0a86899.svg?invert_in_darkmode&sanitize=true" align=middle width=25.95231374999999pt height=29.190975000000005pt/>.

An algorithm for finding <img src="/tex/e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode&sanitize=true" align=middle width=11.027402099999989pt height=22.465723500000017pt/> under the above-mentioned assumptions was proposed by [Peters et al.(2016)]( https://doi.org/10.1111/rssb.12167) for the case when the experimenters already have a sample of the joint distribution of all covariates (i.e. a sample from the *observational distribution*), and samples from different interventions at one or more covariates from <img src="/tex/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode&sanitize=true" align=middle width=14.908688849999992pt height=22.465723500000017pt/>. A concise reimplementation of their method and some examples with synthetic datasets are provided at [[code]](https://github.com/lkania/causal-inference-using-invariant-prediction/blob/master/causality.R)