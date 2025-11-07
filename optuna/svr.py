import optuna
import numpy as np
import pandas as pd
import pickle
import plotly.io as pio

from sklearn.svm import SVR
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
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

# Define objective function
def objective(trial):
    try:
        C = trial.suggest_loguniform("C", 1e-2, 1e3)
        epsilon = trial.suggest_loguniform("epsilon", 1e-3, 1.0)
        kernel = trial.suggest_categorical("kernel", ["linear", "rbf", "poly"])
        gamma = trial.suggest_categorical("gamma", ["scale", "auto"])
        degree = 3 if kernel != "poly" else trial.suggest_int("degree", 2, 5)

        model = SVR(C=C, epsilon=epsilon, kernel=kernel, gamma=gamma, degree=degree)
        #model = make_pipeline(StandardScaler(), svr)

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

# Create multi-objective study
study = optuna.create_study(
    directions=["maximize", "minimize"],
    study_name="svr2_multiobj_log2TPM_150",
    storage="sqlite:///svr2_study_150.db",
    load_if_exists=True
)
study.optimize(objective, n_trials=150)

# Save study
with open("svr2_spearman_mse_study_150.pkl", "wb") as f:
    pickle.dump(study, f)

# Print Pareto front
print("📈 Pareto front trials (SVR):")
for t in study.best_trials:
    print(f"  Corr: {t.values[0]:.4f}, MSE: {t.values[1]:.4f}, Params: {t.params}")

# Plot and save Pareto front
import optuna.visualization as vis
fig = vis.plot_pareto_front(study, target_names=["Spearman", "MSE"])
fig.write_image("optuna_svr2_pareto_front_150.png")
fig.show()

