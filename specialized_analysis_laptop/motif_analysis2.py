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
#regr = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=50)
#regr.fit(X, y)
#print(regr.feature_importances_)
#np.savetxt("/Users/tfriedrich/Downloads/output_regression.txt", regr.feature_importances_)

## names of TFs
#with open('/Users/tfriedrich/Downloads/X_values.xls') as f:
#    first_line = f.readline()


#clf = RandomForestRegressor()
#param_grid = {
#                 'n_estimators': [5, 10, 15, 20],
#                 'max_depth': [2, 5]
#            }
#grid_clf = GridSearchCV(clf, param_grid, cv=10)
#grid_clf.fit(X,y)
#print (grid_clf.best_estimator_)
#model = grid_clf.best_estimator_
#model.fit(X,y_classifier)
model = RandomForestRegressor(bootstrap=True, criterion='mse', max_depth=2,
           max_features='auto', max_leaf_nodes=None,
           min_impurity_decrease=0.0, min_impurity_split=None,
           min_samples_leaf=1, min_samples_split=2,
           min_weight_fraction_leaf=0.0, n_estimators=50, n_jobs=1,
           oob_score=False, random_state=None, verbose=0, warm_start=False)
model.fit(X,y)
np.savetxt("/Users/tfriedrich/Downloads/output_regression5.txt", model.feature_importances_)