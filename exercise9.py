########Question 1########
#load packages and file
import pandas
import scipy
from scipy import optimize
from scipy import stats
from scipy.stats import norm
import numpy
from plotnine import *
#note: used cygwin to change mutation names to numbers
#WT=0, M124K=1, V456D=2, I213N=3
ponzr = pandas.read_csv("ponzr1.csv")
#subset WT and mutation 1
subset1 = ponzr.loc[ponzr.mutation.isin(['0','1']),:]
#subset WT and mutation 2
subset2 = ponzr.loc[ponzr.mutation.isin(['0','2']),:]
#subset WT and mutation 3
subset3 = ponzr.loc[ponzr.mutation.isin(['0','3']),:]

#creating nll function
#null hypothesis
def null(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll= -1*scipy.stats.norm(expected,sigma).logpdf(obs.ponzr1Counts).sum()
    return nll
#alternate hypothesis
def alter(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*obs.mutation
    nll= -1*scipy.stats.norm(expected,sigma).logpdf(obs.ponzr1Counts).sum()
    return nll

#minimizing nll for first treatment (M124K)
initialGuess=numpy.array([2000,300])
alterGuess=numpy.array([2000,-55,300])
fitNull=scipy.optimize.minimize(null,initialGuess,method="Nelder-Mead",options={'disp':True},args=subset1)
fitAlter=scipy.optimize.minimize(alter,alterGuess,method="Nelder-Mead",options={'disp':True},args=subset1)
print(fitNull)
print(fitAlter)
#likelihood ratio test chi2 for M124K
D=2*(fitNull.fun-fitAlter.fun)
#resulting p-value M124K
1-scipy.stats.chi2.cdf(x=D,df=1)


#minimizing nll for V456D
initialGuess=numpy.array([2000,300])
alterGuess=numpy.array([2000,-55,300])
fitNull=scipy.optimize.minimize(null,initialGuess,method="Nelder-Mead",options={'disp':True},args=subset2)
fitAlter=scipy.optimize.minimize(alter,alterGuess,method="Nelder-Mead",options={'disp':True},args=subset2)
print(fitNull)
print(fitAlter)
#likelihood ratio test chi2 for V456D
D=2*(fitNull.fun-fitAlter.fun)
1-scipy.stats.chi2.cdf(x=D,df=1)
#resulting p-value V456D
print(D)

#minimizing nll for I213N
initialGuess=numpy.array([2000,300])
alterGuess=numpy.array([2000,-55,300])
fitNull=scipy.optimize.minimize(null,initialGuess,method="Nelder-Mead",options={'disp':True},args=subset3)
fitAlter=scipy.optimize.minimize(alter,alterGuess,method="Nelder-Mead",options={'disp':True},args=subset3)
print(fitNull)
print(fitAlter)
#likelihood ratio test chi2 I213N
D=2*(fitNull.fun-fitAlter.fun)
1-scipy.stats.chi2.cdf(x=D,df=1)
#p-value for I213N
print(D)

#Question 2#
import pandas
import scipy
