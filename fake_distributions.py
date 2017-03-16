import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns


rand_distribution = np.random.normal(1.0, .5, 10000)


bootstrap_mean = []

for _ in range(200):
    bootstrap_mean.append(np.mean(random.sample(list(rand_distribution), 200)))


fig1, ax = plt.subplots()
#ax.boxplot([rand_distribution,
#            bootstrap_mean])
sns.distplot(rand_distribution)
sns.distplot(bootstrap_mean)
ax.set_yscale('log')
ax.set_xticklabels(["normal", "tumoral"])
ax.set_ylabel("Expression mean distribution")
plt.savefig("bootstrap.png")