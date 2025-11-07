import optuna
import pandas as pd
import numpy as np
import pickle
import plotly.io as pio

from xgboost import XGBRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr

# Load data
traindf = pd.read_csv('/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/control01bigmat_all_rmdupgenes.csv')[[
    'gene_name', 'meancov1k', 'meancovNDR', 'meanFFT950',
    'meanFFTNDR', 'meanSL1k', 'meansimpson1k_100_300', 'TPM'
]]
traindf = traindf.dropna().reset_index(drop=True)

X = traindf.iloc[:, 1:-1]
y = np.log2(traindf['TPM'] + 1)

# Optuna objective
def objective(trial):
    try:
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 100, 1000),
            "max_depth": trial.suggest_int("max_depth", 3, 15),
            "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
            "subsample": trial.suggest_float("subsample", 0.5, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-4, 10.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-4, 10.0, log=True),
            "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
            "verbosity": 0,
            "random_state": 42,
            "n_jobs": -1
        }

        model = XGBRegressor(**params)

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

        return np.mean(spearman_scores), np.mean(mse_scores)

    except Exception as e:
        print(f"Trial failed: {e}")
        return 0.0, float("inf")

# Create and run study
study = optuna.create_study(
    directions=["maximize", "minimize"],
    study_name="xgb_multiobj_log2TPM",
    storage="sqlite:///xgb_study.db",
    load_if_exists=True
)
study.optimize(objective, n_trials=150)

# Save study
with open("xgb_spearman_mse_study.pkl", "wb") as f:
    pickle.dump(study, f)

# Print Pareto front
print("📈 Pareto front trials (XGBoost):")
for t in study.best_trials:
    print(f"  Corr: {t.values[0]:.4f}, MSE: {t.values[1]:.4f}, Params: {t.params}")

# Visualize Pareto front
import optuna.visualization as vis
fig = vis.plot_pareto_front(study, target_names=["Spearman", "MSE"])
fig.write_image("optuna_xgb_pareto_front.png")
fig.show()

