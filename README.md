This package contains the data and the computer code to reproduce the results in 
Zhang, Pan and Wang's paper titled "Community Detection in Attributed Collaboration Network for Statisticians".



It contains 1 data file: 

coauthorship(43journal)_core3.RData     -- the "igraph" graph with 1489 nodes and their attributes



The Code folder contains our own code to reproduce the results in the paper, as well as the R code for ECV (Li et al., 2020), ANCA (Falih et al., 2017), SCORE (Jin et al., 2015; Ji et al., 2016) and the comparison metrics:

main.R     -- the main code to reproduce the results in the paper

ECV.R     -- the functions used in the ECV

ANCA.R    --the functions used in the ANCA

SCORE.R    --the functions used in the SCORE

comparison_metrics.R    --the metrics to measure community detection results





Jin, J. et al. (2015). Fast community detection by score. Annals of Statistics, 43(1):57-89.

Ji, P., Jin, J., et al. (2016). Coauthorship and citation networks for statisticians. The Annals of Applied Statistics, 10(4):1779-1812.

Falih, I., Grozavu, N., Kanawati, R., and Bennani, Y. (2017). ANCA: Attributed network clustering algorithm. In International Conference on Complex Networks and their Applications, pages 241-252. Springer.

Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge sampling. Biometrika, 107(2):257-276.