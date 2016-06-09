from math import exp
def power_option(S,T,sigma,r,n,alpha):
	level=[]
	u=exp(sigma*(T/n)**.5)
	d=1/u
	p=(exp(r*(T/n))-d)/(u-d)
	for i in range(0,n+1):
		pr = S * u**(n-i) * d**i
		vl = pr**alpha
		level.append(vl)
	for i in range(n,0,-1):
		level1=[]
		for j in range(0,i):
			vl = exp(-r*T/n)*(level[j]*p+level[j+1]*(1-p))
			level1.append(vl)
		level = level1
		outputlevel=[]
		if len(level)<=11:
			for i in range(len(level)):
				outputlevel.append(round(level[i],3))
			print (outputlevel)
	for i in range(11,0,-1):
		pr=[]
		for k in range(0,i):
			pr.append(round(S * u**(i-1-k) * d**k,3))
		print (pr)
	return level[0]

value = power_option(100,3/12,0.3,0.1,100,3)
print()
print ("The Power Option price is ",round(value,3))
