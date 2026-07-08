#!/usr/bin/env python3

import pandas as pd
import joblib

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix
)

# ----------------------------------------------------
# Load Data
# ----------------------------------------------------

X = pd.read_csv("selected_features.csv", index_col=0)

labels = pd.read_csv("labels.csv")

y = labels["Label"]

print("Feature Matrix Shape :", X.shape)
print("Labels Shape         :", y.shape)

# ----------------------------------------------------
# Train-Test Split
# ----------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.30,
    random_state=42,
    stratify=y
)

print("\nTraining Samples :", len(X_train))
print("Testing Samples  :", len(X_test))

# ----------------------------------------------------
# Build Random Forest Model
# ----------------------------------------------------

rf = RandomForestClassifier(

    n_estimators=300,
    criterion="gini",
    max_depth=None,
    min_samples_split=2,
    min_samples_leaf=1,
    random_state=42,
    n_jobs=-1

)

# ----------------------------------------------------
# Train
# ----------------------------------------------------

rf.fit(X_train, y_train)

print("\nRandom Forest Model Trained Successfully.")

# ----------------------------------------------------
# Prediction
# ----------------------------------------------------

y_pred = rf.predict(X_test)

# ----------------------------------------------------
# Accuracy
# ----------------------------------------------------

accuracy = accuracy_score(y_test, y_pred)

print("\nAccuracy : {:.2f}%".format(accuracy * 100))

# ----------------------------------------------------
# Classification Report
# ----------------------------------------------------

print("\nClassification Report\n")

print(classification_report(
    y_test,
    y_pred,
    target_names=["ASD","ADHD","IDDD"]
))

# ----------------------------------------------------
# Confusion Matrix
# ----------------------------------------------------

print("\nConfusion Matrix\n")

cm = confusion_matrix(y_test, y_pred)

print(cm)

# ----------------------------------------------------
# Prediction Probabilities
# ----------------------------------------------------

prob = rf.predict_proba(X_test)

prob_df = pd.DataFrame(
    prob,
    columns=["ASD_Prob","ADHD_Prob","IDDD_Prob"]
)

prob_df.to_csv(
    "rf_prediction_probabilities.csv",
    index=False
)

# ----------------------------------------------------
# Feature Importance
# ----------------------------------------------------

importance = pd.DataFrame({

    "Gene": X.columns,

    "Importance": rf.feature_importances_

})

importance = importance.sort_values(
    by="Importance",
    ascending=False
)

importance.to_csv(
    "rf_feature_importance.csv",
    index=False
)

print("\nTop 20 Important Genes\n")

print(importance.head(20))

# ----------------------------------------------------
# Save Predictions
# ----------------------------------------------------

predictions = pd.DataFrame({

    "Actual": y_test.values,

    "Predicted": y_pred

})

predictions.to_csv(
    "rf_predictions.csv",
    index=False
)

# ----------------------------------------------------
# Save Model
# ----------------------------------------------------

joblib.dump(
    rf,
    "multiclass_random_forest.pkl"
)

print("\nModel saved as multiclass_random_forest.pkl")

print("\nPipeline Completed Successfully.")