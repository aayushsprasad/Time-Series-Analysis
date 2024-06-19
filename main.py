import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 5]
from statsmodels.distributions.empirical_distribution import ECDF
import random
from random import seed
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import string
import scipy.stats as stats
#load csv
S = pd.read_csv("SN_d_tot_V2.0.csv")
S.head()
print(S["Daily_Sunspot_Number"].describe())
print(S["Daily_Sunspot_Number"].median())
print(S["Daily_Sunspot_Number"].mode())
print(S["States"].describe())
print(S["States"].median())
print(S["States"].mode())
#Time series plot
series = pd.read_csv('SN_d_tot_V2.0.csv', header=0, index_col=0, parse_dates=True, squeeze=True)
series.plot(title='Daily Number of Sunspots (1818 to 2024)', label='Number of Sunspots')
#Empirical Distribution
y = S['States'].value_counts(normalize=True)
z = y.sort_index(axis=0)
z.plot(title='Fraction of time spent in each State', xlabel='States')
#Plot of States
series = pd.read_csv('States_series.csv', header=0, index_col=0, parse_dates=True, squeeze=True)
series.plot(title='Plot of Markov States', xlabel='Year', ylabel='States')
# Transition Matrix
states = np.array(S['States'])


def TransitionMatrix(states):
    pd.value_counts(states)
    y = len(pd.unique(states))
    P = np.zeros([y, y])

    for i in range(y):
        for j in range(y):
            for x in range(len(states) - 1):
                if states[x] == i and states[x + 1] == j:
                    P[i][j] += 1

    for row in P:
        s = sum(row)
        if s > 0:
            row[:] = [round(f / s, 3) for f in row]
    return P


P = TransitionMatrix(states)
P
#Stationary Distribution
def StationaryDistribution(P):
    A = P.T-np.identity(P.shape[0])
    A = np.vstack([A,np.ones((P.shape[0]))])
    b = np.zeros((P.shape[0])).T
    b = np.zeros((P.shape[0]+1,1))
    b[-1] = 1
    W = np.linalg.lstsq(A,b,rcond=None)[0]
    return W
W = StationaryDistribution(P)
W
#Plot of stationary distribution
series = pd.Series(W.flatten())
series.plot(title='Plot of Stationary Distribution',xlabel='States')
#Markov chain simulation
def SimulateSeries(P,k):
    A = np.zeros(k, dtype=int)
    A[0] = random.choice(range(8))
for t in range(999):
    x = round(random.uniform(0,1),2)
    arr = P_dist[A[t]]
    A[t+1] = min([i for i, e in enumerate(arr) if e >= x])
A
#Time series simulation
def SimulateSeries(P,k):
    P_dist = [np.cumsum(P[i, :]) for i in range(P.shape[0])]
    A = np.zeros(k, dtype=int)
    A[0] = random.choice(range(P.shape[0]))
    for t in range(k-1):
        x = round(random.uniform(0,1),2)
        arr = P_dist[A[t]]
        A[t+1] = min([i for i, e in enumerate(arr) if e >= x])
    return A
A


def TransitionMatrix(states):
    pd.value_counts(states)
    y = len(pd.unique(states))
    P = np.zeros([y, y])

    for i in range(y):
        for j in range(y):
            for x in range(len(states) - 1):
                if states[x] == i and states[x + 1] == j:
                    P[i][j] += 1

    for row in P:
        s = sum(row)
        if s > 0:
            row[:] = [f / s for f in row]
    return P


D = TransitionMatrix(states)
P_dist = [np.cumsum(D[i, :]) for i in range(D.shape[0])]
A = np.zeros(1000, dtype=int)
P_dist
A = np.zeros(1000, dtype=int)
A[0] = random.choice(range(P.shape[0]))
for t in range(999):
    x = round(random.uniform(0,0.95),2)
    print(x)
    arr = P_dist[A[t]]
    print(arr)
    A[t+1] = min([i for i, e in enumerate(arr) if e >= x])
    print([i for i, e in enumerate(arr) if e >= x])
    print(min([i for i, e in enumerate(arr) if e >= x]))

ts1 = pd.Series(states[0:500].flatten(),index=range(500))
ts2 = pd.Series(A[0:500].flatten(), index=range(500))

plt.figure(figsize=(12,5))
plt.xlabel('Time steps')
#plt.ylabel('States')

ax1 = ts1.plot(grid=True, label='Original', lw = 0.7)
ax2 = ts2.plot(color='Orange', grid=True, secondary_y=True, label='Simulated', lw = 0.7)

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()


plt.legend(h1+h2, l1+l2, loc=2)
plt.show()

#Autocorrelation Comparison
def Autocorrelation(X,k):
    X_b = np.average(X)
    n, d = 0, 0
    for i in range(0,len(X) - k):
        n += ((X[i] - X_b)*(X[i+k] - X_b))
    for i in range(0,len(X)):
        d += (X[i] - X_b)**2

    return n/d

states = np.array(S['States'])
ACF1, ACF2 = np.zeros(100), np.zeros(100)
for i in range(len(ACF1)):
    ACF1[i] = Autocorrelation(states,i)
    ACF2[i] = Autocorrelation(A,i)
#import statsmodels
#statsmodels.graphics.tsaplots.plot_acf(ts1,lags=range(500))
df = pd.concat([ts1, ts2], axis=1)
df.corr(method='pearson')

ts1 = pd.Series(ACF1.flatten(),index=range(100))
ts2 = pd.Series(ACF2.flatten(), index=range(100))

plt.figure(figsize=(12,5))
plt.xlabel('k')
#plt.ylabel('ACF')

ax1 = ts1.plot(grid=True, label='Original', lw = 0.7)
ax2 = ts2.plot(color='Orange', grid=False, secondary_y=False, label='Simulated', lw = 0.7)

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()


plt.legend(h1, l1, loc=2)
plt.show()


# Calculating two step Transition Frequencies
def TwoStepFrequency(states):
    pd.value_counts(states)
    y = len(pd.unique(states))
    P = np.zeros([y, y])

    for i in range(y):
        for j in range(y):
            for x in range(len(states) - 2):
                if states[x] == i and states[x + 2] == j:
                    P[i][j] += 1

    for row in P:
        s = sum(row)
        if s > 0:
            row[:] = [round(f / 1, 3) for f in row]
    return P


N = TwoStepFrequency(states)
N

#Two step Transition Frequencies using Transition Matrix
Q = np.round(np.linalg.matrix_power(P, 2),3)
Q

n = np.zeros(N.shape[0])
for i in range(N.shape[0]):
    for j in range(N.shape[0]):
        n[i] += N[i][j]


#GoodnessofFit
n = 0
TS = np.zeros((N.shape[0]))
chi2 = np.zeros(N.shape[0])
p_value = np.zeros(N.shape[0])
for i in range(N.shape[0]):
    N1 = N[i][Q[i]>0]
    n = N1.sum()
    chi2[i] = stats.chi2.ppf(q = 0.95, df = len(N1) - 1)
    for j in range(N.shape[0]):
        if Q[i][j] > 0:
            O = N[i][j]
            E = n * Q[i][j]
            TS[i] += ((O - E)**2 / E)
TS

chi2