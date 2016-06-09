from scipy.stats import norm
from math import log,sqrt,exp
import matplotlib.pyplot as plt
n = norm.pdf
N = norm.cdf
def find_vol(target_value, call_put, S, K, T, r):
    MAX_ITERATIONS = 100
    PRECISION = 10**(-5)
    sigma = 0.5
    for i in range(0, MAX_ITERATIONS):
        price = bs_price(call_put, S, K, T, r, sigma)
        vega = bs_vega(call_put, S, K, T, r, sigma)

        price = price
        diff = target_value - price
        if (abs(diff) < PRECISION):
            return sigma
        sigma = sigma + diff/vega #Xn+1=Xn-f(Xn) / f'(Xn)
    return sigma
def bs_price(cp_flag,S,K,T,r,v,q=0.0):
    d1 = (log(S/K)+(r+v*v/2.)*T)/(v*sqrt(T))
    d2 = d1-v*sqrt(T)
    if cp_flag == 'c':	
        price = S*exp(-q*T)*N(d1)-K*exp(-r*T)*N(d2)
    else:
        price = K*exp(-r*T)*N(-d2)-S*exp(-q*T)*N(-d1)
    return price
def bs_vega(cp_flag,S,K,T,r,v,q=0.0):
    d1 = (log(S/K)+(r+v*v/2.)*T)/(v*sqrt(T))
    return S * sqrt(T)*n(d1)
def bs_calldelta(S,K,T,r,v):
	d1 = (log(S/K)+(r+v*v/2.)*T)/(v*sqrt(T))
	return N(d1)

strike=list(range(91,111))
price=[22.91,22.09,21.27,20.47,19.67,18.88,18.11,17.34,16.59,15.85,15.12,14.41,\
13.71,13.02,12.35,11.69,11.05,10.43,9.82,9.23]
vol=[]
for i in range(len(strike)):
	vol.append(find_vol(price[i],'c',100,strike[i],1,0))
# delta=[]
# for i in range(len(strike)):
	# delta.append(bs_calldelta(100,strike[i],1,0,vol[i]))
# plt.plot(strike,vol)
# plt.ylabel('Implied Volatilities of calls')
# plt.xlabel('Strike')
# plt.show()
# plt.plot(strike,delta)
# plt.ylabel('Delta of calls')
# plt.xlabel('Strike')
# plt.show()
# plt.plot(vol,delta)
# plt.ylabel('Delta of calls')
# plt.xlabel('Implied Volatilities')
# plt.show()

# a=find_vol(bs_price('p',100,100,1,0,0.2,q=0.0), 'p', 100, 90, 1, 0)
# b=find_vol(bs_price('p',100,100,1/12,0,0.2,q=0.0), 'p', 100, 90, 1/12, 0)
# print (a,b)

# print (bs_price('c',100,100,1,0,0.2))



