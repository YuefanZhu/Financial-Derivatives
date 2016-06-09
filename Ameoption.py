# Problem5 Python
# Yuefan Zhu
# yz2835

import math
import matplotlib.pyplot as plt
from scipy.stats import norm
from numpy import log
def price(S,K,T,sigma,r,n,y):
	level=[]
	u=math.exp(sigma*(T/n)**.5)
	d=1/u
	p=(math.exp((r-y)*(T/n))-d)/(u-d)
	for i in range(0,n+1):
		pr = S * u**(n-i) * d**i
		vl = max(pr-K,0)
		level.append(vl)
	for i in range(n,0,-1):
		pracvl=[]
		for k in range(0,i):
			pr = S * u**(i-1-k) * d**k
			vl = max(pr-K,0)
			pracvl.append(vl)
		level1=[]
		for j in range(0,i):
			vl = max(level[j]*p+level[j+1]*(1-p),pracvl[j])/math.exp(r*(T/n))
			level1.append(vl)
		level = level1
		outputlevel=[]
		if len(level)<=11:
			for i in range(len(level)):
				outputlevel.append(round(level[i],2))
			print (outputlevel)
		if len(level)==2:
			delta=(level[0]-level[1])/(S*u-S*d)
	for i in range(11,0,-1):
		pr=[]
		for k in range(0,i):
			pr.append(round(S * u**(i-1-k) * d**k,2))
		print (pr)
	print ("The American Call price is ",round(level[0],2))
	print ("Delta of American Call is ",round(delta,2))
	d1=(log(S*math.exp(-y*T)/(K*math.exp(-r*T)))+0.5*sigma**2*T)/(sigma*T**0.5)
	d2=d1-sigma*T**0.5
	euro=S*math.exp(-y*T)*norm.cdf(d1)-K*math.exp(-r*T)*norm.cdf(d2)
	print ("The European Call price is ",round(euro,2))
	return level

value = price(50,51,1,0.1,0.08,100,0.1)