import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML, display_html
import seaborn as sns
import datetime



from sklearn.model_selection import train_test_split
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import GridSearchCV
from sksurv.nonparametric import kaplan_meier_estimator





data = pd.read_csv('data_2.csv')
data['onset'] = data['onset'].apply(pd.to_datetime, dayfirst=True)
data['confirmation'] = data['confirmation'].apply(pd.to_datetime, dayfirst=True)
data['death_discharge'] = data['death_discharge'].apply(pd.to_datetime, dayfirst=True)
index_censoring = ~data['death_discharge'].isna()
event_censoring_days = []
for i in range(np.shape(data)[0]):
    if index_censoring[i] == False:
        event_censoring_days = np.append(event_censoring_days,
                                         data['confirmation'].dt.dayofyear[i] - data['onset'].dt.dayofyear[i])
    elif index_censoring[i] == True:
        event_censoring_days = np.append(event_censoring_days,
                                         data['death_discharge'].dt.dayofyear[i] - data['onset'].dt.dayofyear[i])
inddaypos = event_censoring_days > 0
event_censoring_days = event_censoring_days[inddaypos]
data = data[inddaypos]
data = data.reset_index()
data = data.drop(['index'], axis=1)


index_censoring = index_censoring[inddaypos]

frm = {'event_censoring_days': event_censoring_days, 'status': index_censoring}
data_y_unstructured = pd.DataFrame(data=frm)
data_y_unstructured = data_y_unstructured[['status', 'event_censoring_days']]
s = data_y_unstructured.dtypes
data_y = np.array([tuple(x) for x in data_y_unstructured.values], dtype=list(zip(s.index, s)))

gender = data['sex']
gender = gender.astype('category')

data.to_csv('data.csv',index=False)
data_y_unstructured.to_csv('data_y.csv',index=False)

median = data['age'].median()


age_type = 'categorized'
age_type = 'median'
if age_type == 'median':
    for i in range(len(data['age'])):
        if data['age'][i] > median:
            data['age'][i] = 'age > 46'
        else:
            data['age'][i] = 'age < 46'
    for value in data['age'].unique():
        mask = data['age'] == value
        time, survival_prob = kaplan_meier_estimator(data_y["status"][mask],data_y["event_censoring_days"][mask])
        survival_prob = survival_prob
        plt.step(time, survival_prob, where="post",
                 label="%s" % (value))
    plt.ylabel("estimated probability of discharge")
    plt.xlabel("time $t$") 
    plt.legend(loc="best")
    limit = plt.gca()
    limit.set_xlim([0,45])
    limit.set_ylim([0, 1])
    plt.savefig('covid_median_age.pdf')
if age_type == 'categorized':
    for i in range(len(data['age'])):
        if data['age'][i] > 60:
            data['age'][i] = 'age > 60'
        elif data['age'][i] <= 60 and data['age'][i] > 46:
            data['age'][i] = '46 < age < 60'
        elif data['age'][i] <= 46 and data['age'][i] > 34:
            data['age'][i] = '34 < age < 46'
        else:
            data['age'][i] = 'age < 34'
    for value in data['age'].unique():
        mask = data['age'] == value
        time, survival_prob =kaplan_meier_estimator(data_y["status"][mask],data_y["event_censoring_days"][mask])
        survival_prob = survival_prob
        plt.step(time, survival_prob, where="post",
                 label="%s" % (value))

    plt.ylabel("Kaplan-Meier estimated probability of discharge")
    plt.xlabel("time $t$")
    plt.legend(loc="best")
    limit = plt.gca()
    limit.set_xlim([0,37])
    limit.set_ylim([0, 1])
    plt.savefig('covid_age_categorized.pdf')

















