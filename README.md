# Rank Regressions for High-dimensional Linear Models with Second-stage Enhancements 
This R package provides functions to estimate the coefficients in high-dimensional linear regressions via a tuning-free and robust approach, Rank Lasso, and second-stage enhancements, Rank SCAD and Rank MCP. To overcome the computational barrier arising from the U-statistics structure in Rank regressions, QICD algorithms and the incomplete U-statistics resampling technique are applied in the package. 

To install the package, please run the following codes in R:

```{r}
library(devtools)
install_github("yunanwu123/RankReg")
```


## Reference

Wang, L., Peng, B., Bradic, J., Li, R. and Wu, Y. (2020), ***A Tuning-free Robust and Efficient Approach to High-dimensional Regression**, Journal of the American Statistical Association, 115:532, 1700-1714*, [doi:10.1080/01621459.2020.1840989](https://doi.org/10.1080/01621459.2020.1840989).

Peng, B. and Wang, L. (2015), ***An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression**, Journal of Computational and Graphical Statistics, 24:3, 676-694*, [doi:10.1080/10618600.2014.913516](https://doi.org/10.1080/10618600.2014.913516).

Clémençon, S., Colin, I., and Bellet, A. (2016), ***Scaling-up empirical risk minimization: optimization of incomplete u-statistics**, The Journal of Machine Learning Research, 17(1):2682–2717*, URL: [https://jmlr.org/papers/v17/15-012.html](https://jmlr.org/papers/v17/15-012.html).

Fan, J. and Li, R. (2001), ***Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties**, Journal of the American Statistical Association, 96:456, 1348-1360*, [doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273). 
