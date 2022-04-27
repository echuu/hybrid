

library('dplyr')                 # data mangaement
options(scipen = 999)
options(dplyr.summarise.inform = FALSE)



# load_libraries.R
# load the packages required to run the algorithm
# library('MASS')                  # misc functions
# library('ggplot2')               # visualization
# library('rpart')                 # rpart()
# library('stringr')               # regex functions
# library('VGAM')                  # for log1mexp() function
# library('matrixcalc')            # trace(), is.positive.definite()
# library('MCMCpack')              # riwish()
# library('CholWishart')           # lmvgamma()
# library('mvtnorm')               # rmvnorm(), dmvnorm()
# library('rpart.plot')            # rpart.rules()
# library('reshape2')              # melt()
# library('rBayesianOptimization') # Matrix_unif()
# library('tmg')
# library('TruncatedNormal')
# library('BayesianTools')
# library('BDgraph')
# library('MLmetrics')


## change value of INSTALL_ALL to TRUE if you are missing the above libraries
INSTALL_ALL = FALSE
if (INSTALL_ALL) {
    install.packages("dplyr")
    install.packages('mvtnorm')
    install.packages('MASS')
    install.packages('ggplot2')
    install.packages('rplot')
    install.packages('tidyr')
    install.packages('readr')
    install.packages('stringr')
    install.packages('dplyr')
    install.packages('VGAM')
    install.packages('matrixcalc')
    install.packages('MCMCpack')
    install.packages('CholWishart')
    install.packages('rpart.plot')
    install.packages('reshape2')
    install.packages('rBayesianOptimization')
    install.packages('tmg')
    install.packages('TruncatedNormal')
    install.packages("BDgraph")
    install.packages("BayesianTools")
}







