# Plackett Luce Random Forest

The <i>plforest</i> function uses the pltree function from the <i>PlackettLuce</i> R package to create a random forest based on bootstrapping and k-fold cross-validation. It takes a ranking class as the response variable and a covariate dataframe as the predictor variable. It predicts the preferred ranking of individuals for a set of objects.

The model setup is as illustrated below:
![plforest](/plforest.png)
