### Imports

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from sklearn.preprocessing import LabelEncoder, StandardScaler, FunctionTransformer
from sklearn.pipeline import make_pipeline

plt.style.use('default')

from sklearn.pipeline import make_pipeline
from sklearn.linear_model import RidgeClassifier
import numpy as np
from fast_borf import BorfBuilder
from fast_borf.pipeline.zero_columns_remover import ZeroColumnsRemover
from fast_borf.pipeline.reshaper import ReshapeTo2D
from fast_borf.pipeline.to_scipy import ToScipySparse
from fast_borf.xai.mapping import BagOfReceptiveFields
from constants import CUSTOM_CONFIG_A3, CUSTOM_CONFIG_A3_NO_DILATION, CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_2_3_4

import logging
logging.basicConfig(filename="log.txt",level=logging.DEBUG)
logging.captureWarnings(True)

### Load data

data = pd.read_csv("temporal_data_with_patient_ready_normalized.csv", index_col=0)
data = data.dropna(axis=0)

columns = data.columns
genes = np.array([name.split("_")[0] for name in list(data.index)])
patients = np.array([name.split("_")[1] for name in list(data.index)])
cells = np.array([name.split("_")[2] for name in list(data.index)])

enc_genes = LabelEncoder()
enc_patients = LabelEncoder()
enc_cells = LabelEncoder()

enc_genes.fit(genes)
enc_patients.fit(patients)
enc_cells.fit(cells)

X = data.values[:, np.newaxis, :]

y_genes = enc_genes.transform(genes)
y_patients = enc_patients.transform(patients)
y_cells = enc_cells.transform(cells)


### Setup the BORF builder
builder = BorfBuilder(
  n_jobs=-2, window_size_min_window_size=None, alphabets_min_symbols=None, alphabets_max_symbols=None, min_window_to_signal_std_ratio=0.15, configs=CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_2_3_4,
  pipeline_objects=[
    (ReshapeTo2D, dict(keep_unraveled_index=True)),
    (ZeroColumnsRemover, dict(axis=0)),
    (ToScipySparse, dict()),
    ],
)
borf = builder.build(X)


### Fit the BORF model
X_transformed = borf.fit_transform(X)
X_transformed_sorted = X_transformed.toarray()[np.argsort(cells)]
X_transformed_sorted = X_transformed_sorted[:, np.argsort(X_transformed.toarray()[np.argsort(cells)].mean(axis=0).ravel())[::-1]]


### Train the classification model
import lightgbm as lgb
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
import xgboost as xgb
from catboost import CatBoostClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.tree import DecisionTreeClassifier
import time


# Define the models to benchmark
models = {
  "LightGBM_OneVsRest": OneVsRestClassifier(lgb.LGBMClassifier(n_jobs=20)),
  "RandomForest_OneVsRest": OneVsRestClassifier(RandomForestClassifier(n_jobs=20)),
  "SVM_OneVsRest": OneVsRestClassifier(SVC(probability=True)),
  "KNN_OneVsRest": OneVsRestClassifier(KNeighborsClassifier(n_jobs=20)),
  "XGBoost_OneVsRest": OneVsRestClassifier(xgb.XGBClassifier(n_jobs=20)),
  "DecisionTree_OneVsRest": OneVsRestClassifier(DecisionTreeClassifier()),
  "RidgeClassifier_OneVsRest": OneVsRestClassifier(RidgeClassifier()),
  "LightGBM": lgb.LGBMClassifier(n_jobs=20),
  "RandomForest": RandomForestClassifier(n_jobs=20),
  "SVM": SVC(probability=True),
  "KNN": KNeighborsClassifier(n_jobs=20),
  "XGBoost": xgb.XGBClassifier(n_jobs=20),
  "DecisionTree": DecisionTreeClassifier(),
  "RidgeClassifier": RidgeClassifier()
}

# Dictionary to store the results
results = {}

# Train and evaluate each model
for model_name, model in models.items():

  start_time = time.time()
  model.fit(X_transformed.astype(float), y_genes)
  y_pred = model.predict(X_transformed.astype(float))
  end_time = time.time()
  
  accuracy = accuracy_score(y_genes, y_pred)
  report = classification_report(y_genes, y_pred, output_dict=True)
  elapsed_time = end_time - start_time
  
  results[model_name] = {
  "accuracy": accuracy,
  "precision": report["weighted avg"]["precision"],
  "recall": report["weighted avg"]["recall"],
  "f1_score": report["weighted avg"]["f1-score"],
  "time": elapsed_time
  }


# Save the results to a CSV file
results_df = pd.DataFrame(results).T
results_df.to_csv("model_benchmark_results.csv")

