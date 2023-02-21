import numpy as np
import matplotlib.pyplot as plt
import random, math

num_sims = 10
mut = 0.01

K1 = 100
K2 = 150
rN = 0.03
rP = 0.06
rM = 0.08
dN = 0.01

tend = 2500
wound = 500

norm_macro = []
pre_macro = []
mal_macro = []
t_macro = []
extinct = []

def cP(t):
    if t<wound:
        return .04
    else:
        return .04*np.tanh(.001*(t-wound))
def cM(t):
    if t<wound:
        return .08
    else:
        return .08*np.tanh(.001*(t-wound))
def sN(t):
    if t>wound and t<wound+2:
        return 1
    else:
        return 0
def sP(t):
    if t>wound and t<wound+2:
        return 1
    else:
        return 0
def sM(t):
    if t>wound and t<wound+2:
        return 1
    else:
        return 0


def run_this():
        normal = [10]
        pre = [0]
        mal = [0]
        t = [0]

        while t[-1] < tend:
            
                current_n = normal[-1]
                current_pre = pre[-1]
                current_mal = mal[-1]

                #broken into birth, death, switching

                rates = [(current_n*rN), current_n*(rN*((current_pre+current_n)/K1)+dN+sN(t[-1])),(current_pre*rP),current_pre*(rP*((current_pre+current_n)/K1)+dN+cP(t[-1])+sP(t[-1])),(current_mal*rM),current_mal*(rM*(current_mal/K2)+dN+cM(t[-1])+sM(t[-1]))]

                rate_sum = sum(rates)

                if rate_sum == 0:
                        extinct.append(t[-1])
                        break

                tau = np.random.exponential(scale=1/rate_sum)

                t.append(t[-1] + tau)

                rand = random.uniform(0,1)

                #Normal cell division event
                if rand * rate_sum <= rates[0]:
                        mal.append(mal[-1])
                        if random.uniform(0,1)>mut:
                            normal.append(normal[-1] + 1)
                            pre.append(pre[-1])
                        else:
                            normal.append(normal[-1])
                            pre.append(pre[-1]+1)


                #Normal cell death event
                elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
                        normal.append(normal[-1] - 1)
                        pre.append(pre[-1])
                        mal.append(mal[-1])

                #Precancerous division event
                elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
                        normal.append(normal[-1])
                        if random.uniform(0,1)>mut:
                            pre.append(pre[-1] + 1)
                            mal.append(mal[-1])
                        else:
                            pre.append(pre[-1])
                            mal.append(mal[-1]+1)
                            
                #Precancerous cell death event
                elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
                        normal.append(normal[-1])
                        pre.append(pre[-1]-1)
                        mal.append(mal[-1])


                #Malignant division event
                elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
                        normal.append(normal[-1])
                        pre.append(pre[-1])
                        mal.append(mal[-1]+1)
                        
                #Malignant cell death event
                elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):
                        normal.append(normal[-1])
                        pre.append(pre[-1])
                        mal.append(mal[-1]-1)       
                            
        norm_macro.append(normal)
        pre_macro.append(pre)
        mal_macro.append(mal)
        t_macro.append(t)

for i in range(num_sims):
        run_this()

lwi = 0.3

extinct_events = len(extinct)
avg_extinct = np.mean(extinct)
std_extinct = np.std(extinct)

print(extinct_events)
print(avg_extinct)
print(std_extinct)

plt.figure()
fig, ax = plt.subplots()
ax.set_title('Acute Wounding')
ax.plot(t_macro[0],norm_macro[0], label='Normal', c='r',lw=lwi)
ax.plot(t_macro[0],pre_macro[0], label='Precancerous', c='b',lw=lwi)
ax.plot(t_macro[0],mal_macro[0], label='Malignant', c='k',lw=lwi)


for i in range(1,num_sims):
    ax.plot(t_macro[i],norm_macro[i], c='r',lw=lwi)
    ax.plot(t_macro[i],pre_macro[i], c='b',lw=lwi)
    ax.plot(t_macro[i],mal_macro[i], c='k',lw=lwi)
    
# get the legend object
leg = ax.legend()

# change the line width for the legend
for line in leg.get_lines():
    line.set_linewidth(2)
    
plt.xlim(0,tend)
plt.ylim(0)
ax = plt.gca()
plt.ylabel('Population Size')
plt.xlabel('Time')
plt.tight_layout()
plt.show()
