# Local Volatility Tree:
from scipy.stats import norm
from math import log,sqrt,exp
import matplotlib.pyplot as plt
n = norm.pdf
N = norm.cdf

def find_vol(target_value, call_put, s, K, T, r):
    MAX_ITERATIONS = 100
    PRECISION = 10**(-5)
    sigma = 0.5
    for i in range(0, MAX_ITERATIONS):
        price = bs_price(call_put, s, K, T, r, sigma)
        vega = bs_vega(call_put, s, K, T, r, sigma)

        price = price
        diff = target_value - price
        if (abs(diff) < PRECISION):
            return sigma
        sigma = sigma + diff/vega #Xn+1=Xn-f(Xn) / f'(Xn)
    return sigma
def bs_price(cp_flag,s,K,T,r,v,q=0.0):
    d1 = (log(s/K)+(r+v*v/2.)*T)/(v*sqrt(T))
    d2 = d1-v*sqrt(T)
    if cp_flag == 'c':	
        price = s*exp(-q*T)*N(d1)-K*exp(-r*T)*N(d2)
    else:
        price = K*exp(-r*T)*N(-d2)-s*exp(-q*T)*N(-d1)
    return price
def bs_vega(cp_flag,s,K,T,r,v,q=0.0):
    d1 = (log(s/K)+(r+v*v/2.)*T)/(v*sqrt(T))
    return s * sqrt(T)*n(d1)



def sigma(S):
	return 0.2*exp(-2*(S/100-1))

def print_tree(tree):
	for i in list(range(len(tree))):
		print(tree[i])

		
		
calls=list()
puts=list()
for K in list(range(95,106)):
		
	price_tree=list()
	prob_tree=list()
	localvol_tree=list()

	T=0.5
	Levels=100
	dt=T/(Levels-1)
	S=100
	price0=list()
	price1=list()
	lvol=list()

	for i in list(range(Levels)):
		
		if (i == 0):
			price1.append(S)
			lvol.append(sigma(S))
			price0=price1
			price_tree.append(price1)
			localvol_tree.append(lvol)
			
		elif ((round((i+1)/2+0.25)-round((i+1)/2+0.25,1))>0):
			price1=list()
			price1.append(S)
			x=S
			for j in list(range(int(i/2))):
				F=price0[int(len(price0)/2)+j]
				sig=lvol[int(len(lvol)/2)+j]
				x=F-F**2*dt*sig**2/(x-F)
				price1.append(x)
			x=S
			for j in list(range(int(i/2))):
				F=price0[int(len(price0)/2)-1-j]
				sig=lvol[int(len(lvol)/2)-1-j]
				x=F+F**2*sig**2*dt/(F-x)
				price1.append(x)
			temp=list()
			for k in list(range(i+1)):
				if (k < i/2):
					temp.append(price1[-1-k])
				else:
					temp.append(price1[k-int(i/2)])
			price1=temp
			price_tree.append(price1)
			lvol=list()
			for k in list(range(i+1)):
				lvol.append(sigma(price1[k]))
			localvol_tree.append(lvol)
			prob=list()
			for k in list(range(i)):
				prob.append((price0[k]-price1[k+1])/(price1[k]-price1[k+1]))
			prob_tree.append(prob)
			price0=price1
		else:
			price1=list()
			price1.append(S*exp(sigma(100)*dt**.5))
			price1.append(S*exp(-sigma(100)*dt**.5))
			x=S*exp(-sigma(100)*dt**.5)
			for j in list(range(int((i-1)/2))):
				F=price0[int((1+len(price0))/2)+j]
				sig=lvol[int((1+len(lvol))/2)+j]
				x=F-F**2*dt*sig**2/(x-F)
				price1.append(x)
			x=S*exp(sigma(100)*dt**.5)
			for j in list(range(int((i-1)/2))):
				F=price0[int((-1+len(price0))/2)-1-j]
				sig=lvol[int((-1+len(lvol))/2)-1-j]
				x=F+F**2*sig**2*dt/(F-x)
				price1.append(x)
			temp=list()
			for k in list(range(i+1)):
				if (k < (i-1)/2):
					temp.append(price1[-1-k])
				else:
					temp.append(price1[k-int((i-1)/2)])
			price1=temp
			price_tree.append(price1)
			lvol=list()
			for k in list(range(i+1)):
				lvol.append(sigma(price1[k]))
			localvol_tree.append(lvol)
			prob=list()
			for k in list(range(i)):
				prob.append((price0[k]-price1[k+1])/(price1[k]-price1[k+1]))
			prob_tree.append(prob)
			price0=price1
	# print('price tree:')
	# print_tree(price_tree)
	# print()
	# print('probability tree:')
	# print_tree(prob_tree)
	# print()
	# print('localvol tree:')
	# print_tree(localvol_tree)

	p=list()
	for i in price_tree[len(price_tree)-1]:
		p.append(max(K-i,0))

	for i in list(range(Levels-1)):
		for j in list(range(Levels-1-i)):
			p[j]=p[j]*prob_tree[Levels-2-i][j]+p[j+1]*(1-prob_tree[Levels-2-i][j])

	c=list()
	for i in price_tree[len(price_tree)-1]:
		c.append(max(i-K,0))

	for i in list(range(Levels-1)):
		for j in list(range(Levels-1-i)):
			c[j]=c[j]*prob_tree[Levels-2-i][j]+c[j+1]*(1-prob_tree[Levels-2-i][j])

	calls.append(c[0])
	puts.append(p[0])
	
print('price tree:')
print_tree(price_tree)
strike=list(range(95,106))
volc=[]
volp=[]
for i in range(len(strike)):
	volc.append(find_vol(calls[i],'c',100,strike[i],T,0))
	volp.append(find_vol(puts[i],'p',100,strike[i],T,0))
plt.plot(strike,volc)
plt.ylabel('Implied Volatilities of calls')
plt.xlabel('Strike')
plt.show()
