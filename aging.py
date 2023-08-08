import numpy as np
import matplotlib.pyplot as plt
import random, math

num_sims = 10
mut = 10**(-3)
mut2 = 10**(-5)

K1 = 10**4
K2 = 10**5
rN = 0.3
rP = 0.6
rM = 0.8
dN = 0.1

tend = 750

norm_macro = []
pre_macro = []
mal_macro = []
t_macro = []

def cP(N,P,M,t):
    eqb = K1-(K1*dN)/rN
    base = 0.3+0.3*np.tanh(P+N+M-eqb)
    return base*np.exp(-.01*t) #set constant for control case

def cM(N,P,M,t):
    eqb = K1-(K1*dN)/rN
    base = 0.6+0.6*np.tanh(P+N+M-eqb)
    return base*np.exp(-.01*t) #set constant for control case

def run_this():
    normal = [int(K1-(K1*dN)/rN)]
    pre = [0]
    mal = [0]
    t = [0]

    while t[-1] < tend:

        current_n = normal[-1]
        current_pre = pre[-1]
        current_mal = mal[-1]

        rates = [current_n * rN * (K1 - current_pre - current_n) / K1, current_n * dN,
                 current_pre * rP * (K1 - current_pre - current_n) / K1, current_pre * (dN + cP(current_n,current_pre,current_mal,t[-1])),
                 current_mal * rM* (K2 - current_mal) / K2, current_mal * (dN + cM(current_n,current_pre,current_mal,t[-1]))]

        rate_sum = sum(rates)

        tau = np.random.exponential(scale=1 / rate_sum)

        t.append(t[-1] + tau)

        rand = random.uniform(0, 1)

        # Normal cell division event
        if rand * rate_sum <= rates[0]:
            mal.append(mal[-1])
            if random.uniform(0, 1) > mut:
                normal.append(normal[-1] + 1)
                pre.append(pre[-1])
            else:
                normal.append(normal[-1])
                pre.append(pre[-1] + 1)


        # Normal cell death event
        elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
            normal.append(normal[-1] - 1)
            pre.append(pre[-1])
            mal.append(mal[-1])

        # Precancerous division event
        elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
            normal.append(normal[-1])
            if random.uniform(0, 1) > mut2:
                pre.append(pre[-1] + 1)
                mal.append(mal[-1])
            else:
                pre.append(pre[-1])
                mal.append(mal[-1] + 1)

        # Precancerous cell death event
        elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
            normal.append(normal[-1])
            pre.append(pre[-1] - 1)
            mal.append(mal[-1])


        # Malignant division event
        elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
            normal.append(normal[-1])
            pre.append(pre[-1])
            mal.append(mal[-1] + 1)

        # Malignant cell death event
        elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):
            normal.append(normal[-1])
            pre.append(pre[-1])
            mal.append(mal[-1] - 1)

    norm_macro.append(normal)
    pre_macro.append(pre)
    mal_macro.append(mal)
    t_macro.append(t)


for i in range(num_sims):
    run_this()
    print(i+1)

lwi = 0.3

fig, ax = plt.subplots()
ax.set_title('Aging')
ax.plot(t_macro[0], norm_macro[0], label='Normal', c='r', lw=lwi)
ax.plot(t_macro[0], pre_macro[0], label='Precancerous', c='b', lw=lwi)
ax.plot(t_macro[0], mal_macro[0], label='Malignant', c='k', lw=lwi)

for i in range(1, num_sims):
    ax.plot(t_macro[i], norm_macro[i], c='r', lw=lwi)
    ax.plot(t_macro[i], pre_macro[i], c='b', lw=lwi)
    ax.plot(t_macro[i], mal_macro[i], c='k', lw=lwi)

leg = ax.legend()
for line in leg.get_lines():
    line.set_linewidth(2)

plt.xlim(0, tend)
plt.ylim(0)
plt.ylabel('Population Size')
plt.xlabel('Time (Days)')
plt.tight_layout()
plt.show()
