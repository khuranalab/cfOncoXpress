import optuna
import numpy as np
import pandas as pd
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, make_scorer
from scipy.stats import spearmanr

# Load and clean data
traindf = pd.read_csv('/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/control01bigmat_all_rmdupgenes.csv')[[
    'gene_name', 'meancov1k', 'meancovNDR', 'meanFFT950',
    'meanFFTNDR', 'meanSL1k', 'meansimpson1k_100_300', 'TPM'
]]
traindf = traindf.dropna().reset_index(drop=True)

X = traindf.iloc[:, 1:-1]
y = np.log2(traindf['TPM'] + 1)

# Define Optuna objective function
def objective(trial):
    try:
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 100, 1000),
            "max_depth": trial.suggest_int("max_depth", 5, 50),
            "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 20),
            "max_features": trial.suggest_categorical("max_features", [None, "sqrt", "log2"]),
            "bootstrap": trial.suggest_categorical("bootstrap", [True, False])
        }

        model = RandomForestRegressor(**params, random_state=42, n_jobs=-1)
        cv = KFold(n_splits=5, shuffle=True, random_state=42)

        spearman_scores = []
        mse_scores = []

        for train_idx, val_idx in cv.split(X):
            X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
            y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

            model.fit(X_train, y_train)
            y_pred = model.predict(X_val)

            spearman_corr, _ = spearmanr(y_val, y_pred)
            mse = mean_squared_error(y_val, y_pred)

            spearman_scores.append(spearman_corr)
            mse_scores.append(mse)

        mean_corr = np.mean(spearman_scores)
        mean_mse = np.mean(mse_scores)

        # Multi-objective return: maximize correlation, minimize MSE
        return mean_corr, mean_mse

    except Exception as e:
        print(f"Trial failed: {e}")
        return 0.0, float("inf")  # Worst-case for both

# Run study — multi-objective
study = optuna.create_study(
    directions=["maximize", "minimize"],
    study_name="rf150_multiobj_log2TPM",
    storage="sqlite:///rf150_multiobj_study.db",
    load_if_exists=True
)
study.optimize(objective, n_trials=150)

# Save the study
with open("rf150_spearman_mse_study.pkl", "wb") as f:
    pickle.dump(study, f)

# Print best trial(s) from Pareto front
print("📈 Pareto front trials:")
for t in study.best_trials:
    print(f"  Corr: {t.values[0]:.4f}, MSE: {t.values[1]:.4f}, Params: {t.params}")

# Save and show visualization
import optuna.visualization as vis
import plotly.io as pio

fig = vis.plot_pareto_front(study, target_names=["Spearman", "MSE"])
fig.write_image("optuna_rf150_pareto_front.png")  # Requires plotly[kaleido]
fig.show()

