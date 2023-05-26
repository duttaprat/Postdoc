import numpy as np
import pandas as pd 
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier




from sklearn.metrics import confusion_matrix, accuracy_score, f1_score, precision_score, recall_score


clf_names = [
    "Nearest_Neighbors_clf",
    "Linear_SVM_clf",
    "RBF_SVM_clf",
    "Gaussian_Process_clf",
    "Decision_Tree_clf",
    "Random_Forest_clf",
    "Neural_Net_clf",
    "AdaBoost_clf",
    "Naive_Bayes_clf",
    "QDA",
]


classifiers = [
    KNeighborsClassifier(n_neighbors = 3),
    SVC(kernel="linear", C=0.025),
    SVC(gamma=2, C=1),
    GaussianProcessClassifier(1.0 * RBF(1.0)),
    DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    MLPClassifier(alpha=1, max_iter=1000),
    AdaBoostClassifier(),
    GaussianNB(),
    QuadraticDiscriminantAnalysis(),
]


## Dataset (Train & Testing) Loading 



data_folder_path = "/home/pdutta/Rekha_LabWork/Enhancer_Prediction_Project/2023/Dataset/80_cut-off/Fixed_length/200/1-1/80-20/Feature_vectors/"
X_train = np.load(data_folder_path+"train_x.npy")
X_test= np.load(data_folder_path+"test_x.npy")
Y_train = np.load(data_folder_path+"train_y.npy")
Y_test =  np.load(data_folder_path+"test_y.npy")


print(X_train.shape , X_test.shape, Y_train.shape , Y_test.shape)



X_train = X_train.astype(float)
Y_train = Y_train.astype(int)
X_test = X_test.astype(float)
Y_test = Y_test.astype(int)



# ## Run the classifiers 

# In[ ]:


for clf_name, classifier in zip(clf_names, classifiers):
    clf = make_pipeline(StandardScaler(), classifier)
    print(clf_name, classifier)
    print("Scalling done")
    clf.fit(X_train, Y_train)
    print("Classifier fitted done")
    Y_pred = clf.predict(X_test)
    cf_matrix = confusion_matrix(Y_test, Y_pred)
    print(cf_matrix)
    group_names = ['True Neg','False Pos','False Neg','True Pos']
    class_accuracy = 100*cf_matrix.diagonal()/cf_matrix.sum(1)
    print(class_accuracy)
    print("Accuracy for class 1: {:.3f}".format(class_accuracy[1])+"\n"             
          "Accuracy for class 0: {:.3f}".format(class_accuracy[0])+"\n"             
          "Overall accuracy: {:.3f}".format(accuracy_score(Y_test, Y_pred)*100)+"\n"             
          "Precision: {:.3f}".format(precision_score(Y_test, Y_pred))+"\n"             
          "Recall: {:.3f}".format(recall_score(Y_test, Y_pred))+"\n"             
          "F-score:  {:.3f}".format(f1_score(Y_test, Y_pred))+"\n")





