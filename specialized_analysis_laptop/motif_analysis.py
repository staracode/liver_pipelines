import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification

#X = pd.read_csv("~/Downloads/X_values.xls", sep="\t")
#y = pd.read_csv("~/Downloads/y_values.xls", sep="\t")

X = np.fromfile("~/Downloads/X_values.xls")
y = np.fromfile("~/Downloads/y_values.xls")

#X, y = make_classification(n_samples=1000, n_features=4, n_informative=2, n_redundant=0, random_state=0, shuffle=False)
print X
print y

clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
clf.fit(X, y)
print(clf.feature_importances_)