import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.grid_search import GridSearchCV
from sklearn.model_selection import train_test_split

y_pandas = pd.read_csv("/Users/tfriedrich/Downloads/y_values.xls", sep="\t")
y_classifier_pandas = pd.read_csv("/Users/tfriedrich/Downloads/y_values_binary.xls", sep="\t")
y= y_pandas['x']
X = np.loadtxt("/Users/tfriedrich/Downloads/X_values.xls", skiprows=1)

print X.shape
print y.shape
y_classifier= y_classifier_pandas['x']

#classification
#X_train, X_test, y_train, y_test = train_test_split(X, y_classifier, test_size=0.33, random_state=42)

##X, y = make_classification(n_samples=1000, n_features=4, n_informative=2, n_redundant=0, random_state=0, shuffle=False)
#clf = RandomForestClassifier(n_estimators=50, max_depth=2, random_state=0)
clf = RandomForestClassifier()
param_grid = {
                 'n_estimators': [5, 10, 15, 20],
                 'max_depth': [2, 5, 7, 9]
            }
grid_clf = GridSearchCV(clf, param_grid, cv=10)
grid_clf.fit(X,y_classifier)
model = grid_clf.best_estimator_
model.fit(X,y_classifier)
#clf_best = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
#            max_depth=7, max_features='auto', max_leaf_nodes=None,
#            min_impurity_decrease=0.0, min_impurity_split=None,
#            min_samples_leaf=1, min_samples_split=2,
#            min_weight_fraction_leaf=0.0, n_estimators=15, n_jobs=1,
#            oob_score=False, random_state=None, verbose=0,
#            warm_start=False)
#clf_best.fit(X, y_classifier)
#np.savetxt("/Users/tfriedrich/Downloads/output_classifier.txt", clf_best.feature_importances_)
np.savetxt("/Users/tfriedrich/Downloads/output_classifier.txt", model.feature_importances_)
#print (grid_clf.best_estimator_)
#RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
#            max_depth=7, max_features='auto', max_leaf_nodes=None,
#            min_impurity_decrease=0.0, min_impurity_split=None,
#            min_samples_leaf=1, min_samples_split=2,
#            min_weight_fraction_leaf=0.0, n_estimators=15, n_jobs=1,
#            oob_score=False, random_state=None, verbose=0,
#            warm_start=False)


#print (grid_clf.best_score_)
#0.794311717861
#print (grid_clf.best_params_)
#{'n_estimators': 15, 'max_depth': 7}
#np.savetxt("/Users/tfriedrich/Downloads/best_estimator.txt", grid_clf.best_estimator_)
#np.savetxt("/Users/tfriedrich/Downloads/best_score.txt", grid_clf.best_score_)
#np.savetxt("/Users/tfriedrich/Downloads/best_params.txt", grid_clf.best_params_)

#regression 
##X, y = make_regression(n_features=4, n_informative=2, random_state=0, shuffle=False)
#regr = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=50)
#regr.fit(X, y)
#print(regr.feature_importances_)
#np.savetxt("/Users/tfriedrich/Downloads/output_regression.txt", regr.feature_importances_)

## names of TFs
#with open('/Users/tfriedrich/Downloads/X_values.xls') as f:
#    first_line = f.readline()

