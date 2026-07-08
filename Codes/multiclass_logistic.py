#!/usr/bin/env python3

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

# Load Data
X = pd.read_csv("selected_features.csv", index_col=0)
labels = pd.read_csv("labels.csv")
y = labels["Label"]
print("Feature Matrix Shape:", X.shape)
print("Number of Classes:", len(set(y)))

# Create Logistic Regression Model
model = LogisticRegression(
    multi_class="multinomial",
    solver="lbfgs",
    max_iter=1000,
    random_state=99
)
model.fit(X, y)

# Predict
predictions = model.predict(X)

# Accuracy
accuracy = accuracy_score(y, predictions)
print("\nAccuracy:", accuracy)

# Confusion Matrix
print("\nConfusion Matrix")
print(confusion_matrix(y, predictions))

# Classification Report
print("\nClassification Report")
print(classification_report(y, predictions))

# Save Predictions
results = pd.DataFrame({
    "Sample": X.index,
    "Actual": y,
    "Predicted": predictions
})
results.to_csv("predictions.csv", index=False)
print("\nPredictions saved.")