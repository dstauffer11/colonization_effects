import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv('data/independence_data.csv', header=0, sep=', ')

fig, ax = plt.subplots(figsize=(9,4))
data.colonizer.value_counts().plot.bar(zorder=5, color='royalblue', ax=ax)
plt.xticks(rotation=0)
ax.yaxis.grid()
ax.set_title('Number of Colonies per Colonizer')
plt.show()


fig, ax = plt.subplots(figsize=(9,4))
data.independence_year.plot.hist(zorder=5, color='royalblue', ax=ax)
ax.yaxis.grid()
ax.set_title('Indepence Year')
plt.show()