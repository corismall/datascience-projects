#intro to machine learning- Kaggle course

'''First familiarize yourself with the data you'll be working with,
save filepath to variable'''

import pandas as pd

melbourne_filepath = '../input/melbourne-housing-snapshot/melb_data.csv'

#read data and store in dataframe

melbourne_data = pd.read_csv(melbourne_filepath)

#print summary of dataframe, count shows how many rows have non-missing values
melbourne_data.describe()

#print list of columns

melbourne_data.columns

#to drop na values

melbourne_data = melbourne_data.dropna(axis=0)

###############################################
'''To subset your dataframe: use dot notation to select prediction target ("y")
You can pull out a variable with dot-notation. This single column is stored in
a Series, which is broadly like a DataFrame with only a single column of data.'''

y = melbourne_data.Price

#selecting features with a column list ("X")

melbourne_features = ['list of fields']
X = melbourne_data[melbourne_features]

#########################################################
''' scikit-learn library is used to create model (sklearn), the most popular library
for modeling data in dataframe structure. Here are the steps for creating your
model'''

#Step1: define model
#Step2: fit model to data (heart of modeling)
#Step3: predict
#Step4: evaluate output

from sklearn.tree import DecisionTreeRegressor

#Step1: specify a number for random_state to ensure same results each run

melbourne_model = DecisionTreeRegressor(random_state=1)

#Step2: fit your data

melbourne_model.fit(X, y)

print("Making predictions for the following 5 houses:")
print(X.head())
print("The predictions are")
print(melbourne_model.predict(X.head()))

##############################################

'''How to measure model quality "Model Validation":
There are many metrics for summarizing model quality, but we'll start with
one called Mean Absolute Error (also called MAE) MAE = mean(abs(actual−predicted))
on average, our predictions are off by about MAE '''

from sklearn.metrics import mean_absolute_error

predicted_home_prices = melbourne_model.predict(X)
mean_absolute_error(y, predicted_home_prices)

'''using a single "sample" of items for both building the model and evalutating
it generates "in sample" scores, "in sample" scores give the false impression
that the model is accurate --- its important to evaluate the model using new data
termed "Validation Data" to get an "out of sample" score

comparing the "in sample" to the "out of sample" scores will give you the difference
between an almost exactly right model and one that is likely practically unusable'''

'''a simple way to generate Validation data is to split the training data so the model
doesn't see all available data at first'''

from sklearn.model_selection import train_test_split

'''split data into training and validation data, for both features and target
The split is based on a random number generator. Supplying a numeric value to
the random_state argument guarantees we get the same split every time you
run this script.'''

train_X, val_X, train_y, val_y = train_test_split(X, y, random_state = 0)

melbourne_model = DecisionTreeRegressor()
melbourne_model.fit(train_X, train_y)

#predict prices using validation data
val_predictions = melbourne_model.predict(val_X)
print(mean_absolute_error(val_y, val_predictions))

#################################################
