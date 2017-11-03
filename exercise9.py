########Question 1########
#load packages and file
import pandas
import scipy
from scipy import optimize
from scipy import stats
from scipy.stats import norm
import numpy
from plotnine import *
#note: used sed in cygwin to change mutation names to numbers
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
import numpy
import scipy
from scipy.optimize import minimize
from scipy.stats import norm
mgrowth = open("MmarinumGrowth.csv", "r")
def (y,S,umax,Ks):
    sigma = y[0]
    u = umax*S/(S+Ks)
    nll = -1*norm(u,sigma).logpdf(S.u).sum()
    return nll
guess=numpy.array([1,1,1])
fit=minimize(sim,guess,method='Nelder-Mead',options={'disp': True},args=17)
#####Question 3#####
#load packages
import pandas
import scipy
from scipy import optimize
from scipy import stats
from scipy.stats import norm
import numpy
from plotnine import *
#load dataset
leafy = pandas.read_csv("leafDecomp.csv")

#define custom likelihood functions
def constant(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll= -1*scipy.stats.norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll
    
#alternate hypotheses
def linear(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*obs.Ms
    nll= -1*scipy.stats.norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll
    
def humped(p,obs):
    B0=p[0]
    B1=p[1]
    B2=p[2]
    sigma=p[3]
    
    expected=B0+B1*obs.Ms+B2*(obs.Ms)**2
    nll= -1*scipy.stats.norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

#parameter estimates
###NOTE TO MATI- HUMPED FUNCTION DOESN'T WORK CURRENTLY####
constantGuess=numpy.array([1,1])
linearGuess=numpy.array([1,1,1])
humpedGuess=numpy.array([200,10,-0.2,1])
fitConstant=scipy.optimize.minimize(constant,constantGuess,method="Nelder-Mead",options={'disp':True},args=leafy)
fitLinear=scipy.optimize.minimize(linear,linearGuess,method="Nelder-Mead",options={'disp':True},args=leafy)
fitHumped=scipy.optimize.minimize(humped,humpedGuess,method="Nelder-Mead",options={'disp':True},args=leafy)




