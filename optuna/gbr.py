import optuna
import pandas as pd
import numpy as np
import pickle
import plotly.io as pio

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr

# Load and clean data
traindf = pd.read_csv('/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/control01bigmat_all_rmdupgenes.csv')[[
    'gene_name', 'meancov1k', 'meancovNDR', 'meanFFT950',
    'meanFFTNDR', 'meanSL1k', 'meansimpson1k_100_300', 'TPM'
]]
traindf = traindf.dropna().reset_index(drop=True)

X = traindf.iloc[:, 1:-1]
y = np.log2(traindf['TPM'] + 1)

# Optuna objective for GradientBoostingRegressor
def objective(trial):
    try:
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 100, 1000),
            "max_depth": trial.suggest_int("max_depth", 3, 20),
            "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 20),
            "learning_rate": trial.suggest_loguniform("learning_rate", 1e-3, 0.3),
            "subsample": trial.suggest_float("subsample", 0.5, 1.0),
            "max_features": trial.suggest_categorical("max_features", ["auto", "sqrt", "log2", None])
        }

        model = GradientBoostingRegressor(**params, random_state=42)

        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        spearman_scores = []
        mse_scores = []

        for train_idx, val_idx in cv.split(X):
            X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
            y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

            model.fit(X_train, y_train)
            preds = model.predict(X_val)

            spearman_scores.append(spearmanr(y_val, preds).correlation)
            mse_scores.append(mean_squared_error(y_val, preds))

        mean_corr = np.mean(spearman_scores)
        mean_mse = np.mean(mse_scores)

        return mean_corr, mean_mse

    except Exception as e:
        print(f"Trial failed: {e}")
        return 0.0, float("inf")

# Create study
study = optuna.create_study(
    directions=["maximize", "minimize"],
    study_name="gbr150_multiobj_log2TPM",
    storage="sqlite:///gbr150_study.db",
    load_if_exists=True
)
study.optimize(objective, n_trials=150)

# Save study to file
with open("gbr150_spearman_mse_study.pkl", "wb") as f:
    pickle.dump(study, f)

# Print Pareto front trials
print("📈 Pareto front trials (GBR):")
for t in study.best_trials:
    print(f"  Corr: {t.values[0]:.4f}, MSE: {t.values[1]:.4f}, Params: {t.params}")

# Plot and save figure
import optuna.visualization as vis
fig = vis.plot_pareto_front(study, target_names=["Spearman", "MSE"])
fig.write_image("optuna_gbr150_pareto_front.png")
fig.show()

