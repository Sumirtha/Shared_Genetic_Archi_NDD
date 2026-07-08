#!/usr/bin/env python3

import pandas as pd

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.feature_selection import mutual_info_classif

# ---------------------------------
# Read Feature Matrix
# ---------------------------------

X = pd.read_csv("/home/ibab/PyCharmMiscProject/mrop/feature_matrix.csv", index_col=0)

labels = pd.read_csv("/home/ibab/PyCharmMiscProject/mrop/labels.csv")

y = labels["Label"]

print("Feature Matrix Shape :", X.shape)

# ---------------------------------
# Remove genes with no mutations
# ---------------------------------

X = X.loc[:, (X.sum(axis=0) > 0)]

print("After removing empty genes :", X.shape)

# ---------------------------------
# CHI-SQUARE FEATURE SELECTION
# ---------------------------------

k = 100        # Number of genes to keep

chi_selector = SelectKBest(score_func=chi2, k=min(k, X.shape[1]))

X_chi = chi_selector.fit_transform(X, y)

selected_genes = X.columns[chi_selector.get_support()]

chi_scores = chi_selector.scores_

chi_df = pd.DataFrame({
    "Gene": X.columns,
    "ChiScore": chi_scores
})

chi_df = chi_df.sort_values(
    by="ChiScore",
    ascending=False
)

chi_df.to_csv("chi_square_scores.csv", index=False )

print("Chi-square scores saved.")

# ---------------------------------
# MUTUAL INFORMATION
# ---------------------------------

# mi_scores = mutual_info_classif( X, y, random_state=99 )
#
# mi_df = pd.DataFrame({
#
#     "Gene": X.columns,
#
#     "MI_Score": mi_scores
#
# })
#
# mi_df = mi_df.sort_values( by="MI_Score", ascending=False )
#
# mi_df.to_csv("/mutual_information_scores.csv", index=False )
#
# print("Mutual Information scores saved.")

# ---------------------------------
# Save Selected Features
# ---------------------------------

X_selected = X[selected_genes]

X_selected.to_csv( "selected_features.csv" )

print("Selected feature matrix saved.")

print("\nTop Selected Genes")

print(selected_genes)