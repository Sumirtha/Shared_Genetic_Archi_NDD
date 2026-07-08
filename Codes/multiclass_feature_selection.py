#!/usr/bin/env python3

import pandas as pd

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

# --------------------------------------------
# Load Feature Matrix
# --------------------------------------------

X = pd.read_csv("feature_matrix.csv", index_col=0)

labels = pd.read_csv("labels.csv")

y = labels["Label"]

print("Feature Matrix Shape :", X.shape)

# --------------------------------------------
# Remove genes with no variation
# --------------------------------------------

X = X.loc[:, X.sum(axis=0) > 0]

print("After removing empty genes :", X.shape)

# --------------------------------------------
# Chi-square Feature Selection
# --------------------------------------------

k = min(200, X.shape[1])

selector = SelectKBest(score_func=chi2, k=k)

X_selected = selector.fit_transform(X, y)

selected_columns = X.columns[selector.get_support()]

scores = selector.scores_

# --------------------------------------------
# Save Chi-square Scores
# --------------------------------------------

score_table = pd.DataFrame({

    "Gene": X.columns,
    "Chi2_Score": scores

})

score_table = score_table.sort_values(
    by="Chi2_Score",
    ascending=False
)

score_table.to_csv(
    "chi_square_scores.csv",
    index=False
)

print("Chi-square scores saved.")

# --------------------------------------------
# Selected Feature Matrix
# --------------------------------------------

selected_df = pd.DataFrame(
    X_selected,
    columns=selected_columns,
    index=X.index
)

selected_df.to_csv("selected_features.csv")

print("Selected feature matrix saved.")

# --------------------------------------------
# Save Selected Gene List
# --------------------------------------------

selected_gene_table = pd.DataFrame({
    "Gene": selected_columns
})

selected_gene_table.to_csv(
    "selected_genes.csv",
    index=False
)

# --------------------------------------------
# Display Results
# --------------------------------------------

print("\nTop 20 Selected Genes\n")

print(score_table.head(20))

print("\nSelected Feature Matrix Shape :")

print(selected_df.shape)