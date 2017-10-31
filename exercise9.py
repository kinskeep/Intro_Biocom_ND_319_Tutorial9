#Question 1
#load packages and file
import pandas
import scipy
from scipy import optimize
from scipy import stats
from scipy.stats import norm
import numpy
ponzr = pandas.read_csv("ponzr1.csv")
#subset WT and mutation 1
subset1 = ponzr.loc[ponzr.mutation.isin(['WT','M124K']),:]
#subset WT and mutation 2
subset2 = ponzr.loc[ponzr.mutation.isin(['WT','V456D']),:]
#subset WT and mutation 3
subset3 = ponzr.loc[ponzr.mutation.isin(['WT','I213N']),:]
#creating nll
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*ponzr.ponzr1Counts
    nll= -1*scipy.stats.norm(expected,sigma).logpdf(ponzr.mutation).sum()
    return nll
    print(nllike)
    
#minimizing nll
initialGuess=numpy.array([1,1,1])
fitNull=scipy.optimize.minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp':True},args=1)
fitAlter=scipy.optimize.minimize(nllike,intialGuess,method="Nelder-Mead",options={'disp':True},args=q1)

#likelihood ratio test chi2

D=2*(fitNull.fun-fitAlter.fun)
1-scipy.stats.chi2.cdf(x=D,df=1)


#Question 2#
import pandas
import scipy
