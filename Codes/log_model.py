#!/usr/bin/env python3

import pandas as pd
import joblib

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix
)

# --------------------------------------------------
# Load Data
# --------------------------------------------------

X = pd.read_csv("selected_features.csv", index_col=0)

labels = pd.read_csv("labels.csv")

y = labels["Label"]

print("Feature Matrix Shape :", X.shape)
print("Labels Shape         :", y.shape)

# --------------------------------------------------
# Train-Test Split
# --------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.30,
    random_state=42,
    stratify=y
)

print("\nTraining Samples :", X_train.shape[0])
print("Testing Samples  :", X_test.shape[0])

# --------------------------------------------------
# Build Multiclass Logistic Regression Model
# --------------------------------------------------

model = LogisticRegression(
    multi_class="multinomial",
    solver="lbfgs",
    max_iter=5000,
    random_state=42
)

# --------------------------------------------------
# Train Model
# --------------------------------------------------

model.fit(X_train, y_train)

print("\nModel trained successfully.")

# --------------------------------------------------
# Prediction
# --------------------------------------------------

y_pred = model.predict(X_test)

# --------------------------------------------------
# Accuracy
# --------------------------------------------------

accuracy = accuracy_score(y_test, y_pred)

print("\nAccuracy : {:.2f}%".format(accuracy * 100))

# --------------------------------------------------
# Classification Report
# --------------------------------------------------

print("\nClassification Report\n")

print(classification_report(
    y_test,
    y_pred,
    target_names=["ASD", "ADHD", "IDDD"]
))

# --------------------------------------------------
# Confusion Matrix
# --------------------------------------------------

print("\nConfusion Matrix\n")

cm = confusion_matrix(y_test, y_pred)

print(cm)

# --------------------------------------------------
# Predicted Probabilities
# --------------------------------------------------

probabilities = model.predict_proba(X_test)

prob_df = pd.DataFrame(
    probabilities,
    columns=["ASD_Prob", "ADHD_Prob", "IDDD_Prob"]
)

prob_df.to_csv("prediction_probabilities.csv", index=False)

print("\nPrediction probabilities saved.")

# --------------------------------------------------
# Save Predictions
# --------------------------------------------------

prediction_df = pd.DataFrame({
    "Actual": y_test.values,
    "Predicted": y_pred
})

prediction_df.to_csv("predictions.csv", index=False)

print("Predictions saved.")

# --------------------------------------------------
# Save Model
# --------------------------------------------------

joblib.dump(model, "multiclass_logistic_regression.pkl")

print("Model saved as multiclass_logistic_regression.pkl")

# --------------------------------------------------
# Logistic Regression Coefficients
# --------------------------------------------------

coef_df = pd.DataFrame(
    model.coef_,
    columns=X.columns,
    index=["ASD", "ADHD", "IDDD"]
)

coef_df.to_csv("logistic_regression_coefficients.csv")

print("Coefficient table saved.")

print("\nPipeline Completed Successfully.")