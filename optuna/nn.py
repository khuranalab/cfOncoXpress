import optuna
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# Load data
traindf = pd.read_csv('/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/control01bigmat_all_rmdupgenes.csv')[[
    'gene_name', 'meancov1k', 'meancovNDR', 'meanFFT950',
    'meanFFTNDR', 'meanSL1k', 'meansimpson1k_100_300', 'TPM'
]]
traindf = traindf.dropna().reset_index(drop=True)

X = traindf.iloc[:, 1:-1].values
y = np.log2(traindf['TPM'] + 1).values

# Train/test split
X_train, X_val, y_train, y_val = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Standardize features (needed for NN)
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val = scaler.transform(X_val)

def objective(trial):
    tf.keras.backend.clear_session()  # clear previous models from memory

    model = keras.Sequential()
    n_layers = trial.suggest_int("n_layers", 1, 3)
    for i in range(n_layers):
        units = trial.suggest_int(f"units_l{i}", 16, 128, step=16)
        activation = trial.suggest_categorical(f"act_l{i}", ["relu", "tanh"])
        model.add(layers.Dense(units, activation=activation))
        dropout = trial.suggest_float(f"dropout_l{i}", 0.0, 0.5)
        if dropout > 0:
            model.add(layers.Dropout(dropout))

    model.add(layers.Dense(1))  # Output layer

    lr = trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True)
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=lr),
        loss="mse",
        metrics=["mae"]
    )

    batch_size = trial.suggest_categorical("batch_size", [16, 32, 64])
    history = model.fit(
        X_train, y_train,
        validation_data=(X_val, y_val),
        epochs=50,
        batch_size=batch_size,
        verbose=0,
        callbacks=[keras.callbacks.EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)]
    )

    val_pred = model.predict(X_val).flatten()

    # Calculate Spearman correlation and MSE for multi-objective
    spearman_corr = spearmanr(y_val, val_pred).correlation
    mse = mean_squared_error(y_val, val_pred)

    # Optuna expects tuple for multi-objective (maximize corr, minimize mse)
    return spearman_corr, mse

# Create multi-objective study: maximize Spearman, minimize MSE
study = optuna.create_study(
    directions=["maximize", "minimize"],
    study_name='nn150_multiobj_log2TPM',
    storage="sqlite:///optuna_nn150_study.db",
    load_if_exists=True
)

study.optimize(objective, n_trials=150)

# Print Pareto front trials
print("📈 Pareto front trials (NN):")
for t in study.best_trials:
    print(f"Spearman: {t.values[0]:.4f}, MSE: {t.values[1]:.4f}, Params: {t.params}")

# Save study for later
import pickle
with open("nn150_multiobj_study.pkl", "wb") as f:
    pickle.dump(study, f)

# Visualization: Optimization history of Spearman and MSE
# import optuna.visualization as vis

# fig = vis.plot_optimization_history(study)
# fig.write_image("optuna_optimization_history_nn_150.png")
# fig.show()

# fig_slice = vis.plot_slice(study)
# fig_slice.write_image("optuna_slice_nn_150.png")
# fig_slice.show()

import optuna.visualization as vis
fig = vis.plot_pareto_front(study, target_names=["Spearman", "MSE"])
fig.write_image("optuna_nn150_pareto_front.png")
fig.show()

