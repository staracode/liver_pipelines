import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.grid_search import GridSearchCV
from sklearn.model_selection import train_test_split

y_pandas = pd.read_csv("/Users/tfriedrich/Downloads/y_values.xls", sep="\t")
y= y_pandas['x']
X = np.loadtxt("/Users/tfriedrich/Downloads/X_values.xls", skiprows=1)

print X.shape
print y.shape



#regression 
regr = RandomForestRegressor(max_depth=3, random_state=0, n_estimators=10)
regr.fit(X, y)
print(regr.feature_importances_)
np.savetxt("/Users/tfriedrich/Downloads/output_regression_test2.txt", regr.feature_importances_)
