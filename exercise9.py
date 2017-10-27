#Question 1
import pandas
import scipy
ponzr = pandas.read_csv("ponzr1.csv")
#subset WT and mutation 1
subset1 = ponzr.loc[ponzr.mutation.isin(['WT','M124K']),:]
#subset WT and mutation 2
subset2 = ponzr.loc[ponzr.mutation.isin(['WT','V456D']),:]
#subset WT and mutation 3
subset3 = ponzr.loc[ponzr.mutation.isin(['WT','I213N']),:]
initialGuess=numpy.array([1,1,1])
fitNull=minimize(fun1a,initialGuess,method="Nelder-Mead",options={'disp':True},args=q1)
fitAlter=minimize(fun1b,intialGuess,method="Nelder-Mead",options={'disp':True},args=q1)
1-scipy.stats.chi2.cdf(x=D,df=1)


#Question 2#
import pandas
import scipy
