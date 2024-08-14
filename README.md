# SPARforest
Code for SPARforest

# Paper & Citation

MacBride, C., Davies, V., & Lee, D. (2023). [Conditional autoregressive models fused with random forests to improve small-area spatial prediction](https://arxiv.org/abs/2312.12106). _arXiv preprint arXiv:2312.12106_.

```
@article{macbride2023conditional,
  title={Conditional autoregressive models fused with random forests to improve small-area spatial prediction},
  author={MacBride, Cara and Davies, Vinny and Lee, Duncan},
  journal={arXiv preprint arXiv:2312.12106},
  year={2023}
}
```

# Instructions

The functions for CAR-Forest and SAR-Forest are contained in the R files _CAR-Forest.R_ and _SAR-Forest.R_. Consider a training set with _K_ observations and _p_ features, and a test set with _Q=N-K_ observations and _p_ features. Then the arguments to the _CARForest()_ / _SARForest()_ functions are as follows.

- _X.training_ - The _K x p_ design matrix of features for the training set areal units.
- _X.test_ - The _Q x p_ design matrix of features for the test set areal units. This should contain the same features as _X.training_ but be for the _Q_ test set areal units. 
- _Y.training_ - The _K x 1_ vector of target variables for the training set areal units.
- _coords.training_ - The _K x 2_ matrix of spatial coordinates for the training set areal units (centroids).
- _coords.test_ - The _Q x 2_ matrix of spatial coordinates for the test set areal units (centroids).
- _n_tr_ - The number of trees to use when fitting the random forest component.
- _m_try_ - The value of the $m_{try}$ parameter to use when fitting the random forest component.
- _min_node_ - The value of the $min_{node}$ parameter to use when fitting the random forest component.
- _D_ - The number of nearest neighbours to use when constructing the binary neighbourhood matrix used in the CAR / SAR model component.
- _R_ - The number of iterations of the SPAR-Forest algorithm.
- _logscale_ - Logical - Is the model fitted on the logscale and hence the predictions and 95\% prediction intervals need to be exponentiated before being returned? 
- _inla_samples_ - The number of samples to draw within the INLA software from the CAR / SAR model's posterior predictive distribution for each test set data point. 

Once the model has run it returns a _list_ object with the following 2 elements. 

- _predictions_ - This is a _list_ object with $R$ elements. The $r$th element of this list (for $r=1,\ldots,R$) contains the results from running the algorithm $r$ times, i.e., setting $R=r$. Specifically, it contains a _data.frame_ object with $Q$ rows and 3 columns, where each row relates to an areal unit in the test set in the same order as the _X.test_ and _coords.test_ arguements. The first column contains the predictions, while the second and third columns contain the lower and upper ends of the 95\% prediction intervals. 
- _components_ - This is also a _list_ object with $R$ elements, where again the  _r_ th element of this list  contains the results from running the algorithm $r$ times. Specifically, the _r_ th element of the list contains a _data.frame_ object with $N$ rows and 2 columns, where each row relates to an areal unit which are ordered as the training set observations followed by the test set observations. The two columns relate to the _r_ th estimates of the random forest component $((\hat{m}^{(k)} [x(A_k)])^N_{k=1})$ and the random effect component $((\hat\phi (A_k))_{k=1}^{N})$

Thus, the algorithm only needs to be run once with a given _R_ value to get the results for all _R_ values less than or equal to the one chosen. 
