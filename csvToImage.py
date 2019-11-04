import csv
import matplotlib.pyplot as plt
import numpy as np
modelName = input("モデル名:")
with open('./QA/exp/'+modelName+'.csv')as f:
    reader=csv.reader(f)
    l = [row for row in reader]
    x = [row[0] for row in l]
    y = [row[1] for row in l]

plt.scatter(x,y, label=modelName)
plt.show()
