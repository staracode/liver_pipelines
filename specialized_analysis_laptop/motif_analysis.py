import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression

#X = pd.read_csv("/Users/tfriedrich/Downloads/X_values.xls", sep="\t")
y_pandas = pd.read_csv("/Users/tfriedrich/Downloads/y_values.xls", sep="\t")
y_classifier_pandas = pd.read_csv("/Users/tfriedrich/Downloads/y_values_binary.xls", sep="\t")

X = (np.loadtxt("/Users/tfriedrich/Downloads/X_values.xls", skiprows=1))
#y = np.fromfile("/Users/tfriedrich/Downloads/y_values.xls")
y= y_pandas['x']
#print X.shape
print y.shape
y_classifier= y_classifier_pandas['x']

#X, y = make_classification(n_samples=1000, n_features=4, n_informative=2, n_redundant=0, random_state=0, shuffle=False)
clf = RandomForestClassifier(n_estimators=50, max_depth=2, random_state=0)
clf.fit(X, y_classifier)

#X, y = make_regression(n_features=4, n_informative=2, random_state=0, shuffle=False)
regr = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=50)
regr.fit(X, y)
print(regr.feature_importances_)

#with open('/Users/tfriedrich/Downloads/X_values.xls') as f:
#    first_line = f.readline()

np.savetxt("/Users/tfriedrich/Downloads/output_regression.txt", regr.feature_importances_)
np.savetxt("/Users/tfriedrich/Downloads/output_classifier.txt", clf.feature_importances_)