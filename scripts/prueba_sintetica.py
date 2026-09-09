"""
prueba con datos falsos para ver si el pipeline de pam50 esta bien montado
genera datos falsos con señal conocida, les pasa el mismo pipeline que la
memoria (log1p, top 3000 genes, z-score, random forest) y mira si el modelo
la recupera o no sin tocar la base de datos
"""

import numpy as np
from collections import Counter

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.metrics import accuracy_score, f1_score, classification_report

RANDOM_STATE = 42
N_SAMPLES = 221
N_GENES = 19_938
TOP_GENES = 3000
SUBTIPOS = ["LumA", "LumB", "Basal", "Normal", "Her2"]
N_MARCADORES_POR_CLASE = 200  # num. genes falsos marcan cada subtipo
SIG_LOG2FC = 2.0 # cuanto sube la senial en su clase
RUIDO_BASE = 0.8 # cuanta dispersion de fondo metemos

rng = np.random.default_rng(RANDOM_STATE)


# etiquetas, mas o menos como en la vida real
labels = np.array(
    ["LumA"] * 106 + ["LumB"] * 39 + ["Basal"] * 31 + ["Normal"] * 30 + ["Her2"] * 15
)
print(f"Muestras sintéticas: {len(labels)}")
print(f"Distribución: {dict(sorted(Counter(labels).items()))}")

# matriz de expresion falsa (tpm, todo positivo)
# valores en log con media 3 y un poco de ruido
log_base = rng.normal(loc=3.0, scale=RUIDO_BASE, size=(N_SAMPLES, N_GENES))

# a cada subtipo le subimos sus genes y se los bajamos a los demas
marcador_idx = np.arange(len(SUBTIPOS) * N_MARCADORES_POR_CLASE)
for c, subtipo in enumerate(SUBTIPOS):
    filas_c = np.where(labels == subtipo)[0]
    filas_resto = np.where(labels != subtipo)[0]
    genes_c = marcador_idx[c * N_MARCADORES_POR_CLASE:(c + 1) * N_MARCADORES_POR_CLASE]
    log_base[np.ix_(filas_c, genes_c)] += SIG_LOG2FC # subimos su clase
    log_base[np.ix_(filas_resto, genes_c)] -= SIG_LOG2FC * 0.5  # bajamos al resto

tpm = np.exp(log_base)  # pasamos a tpm sin negativos
print(f"Matriz TPM sintética: {tpm.shape[0]} muestras x {tpm.shape[1]} genes")


# mismo pipeline que en la memoria
X_log = np.log1p(tpm.astype(np.float32))

# top 3000 genes por varianza (antes del z-score)
varianzas = X_log.var(axis=0)
top_idx = np.argsort(varianzas)[::-1][:TOP_GENES]
X_sel = X_log[:, top_idx]

# z-score por gen
media = X_sel.mean(axis=0)
std = X_sel.std(axis=0)
std[std == 0] = 1.0
X_z = (X_sel - media) / std
print(f"Genes seleccionados por varianza: {len(top_idx)} z-score aplicado")

# entrenar y evaluar igual que en la memoria
rf = RandomForestClassifier(
    n_estimators=300, class_weight="balanced", n_jobs=-1, random_state=RANDOM_STATE
)

# split 80/20 estratificado
X_tr, X_te, y_tr, y_te = train_test_split(
    X_z, labels, test_size=0.2, stratify=labels, random_state=RANDOM_STATE
)
rf.fit(X_tr, y_tr)
y_pred = rf.predict(X_te)
acc = accuracy_score(y_te, y_pred)
f1 = f1_score(y_te, y_pred, average="macro")
print("\nTest interno 20% (sintético):")
print(f"Acc.={acc:.4f} F1 macro={f1:.4f}")
print(classification_report(y_te, y_pred, zero_division=0))

# validacion cruzada 5-fold
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
cv_scores = cross_val_score(rf, X_z, labels, cv=cv, scoring="f1_macro", n_jobs=-1)
print(f"CV F1 macro 5-fold: {cv_scores.mean():.4f} +- {cv_scores.std():.4f}")
print(f"Scores por fold: {np.round(cv_scores, 4)}")


# resutados
azar = 1 / len(SUBTIPOS)

print("\nResultado pruebas sinteticas pipeline pam50:\n")
print(f"Nivel de azar (5 clases): F1 ≈ {azar:.2f}")
print(f"Test interno Acc.={acc:.4f}  F1 macro={f1:.4f}")
print(f"CV 5-fold F1 macro={cv_scores.mean():.4f} ± {cv_scores.std():.4f}\n")

if cv_scores.mean() > azar + 0.3:
    print("OK! el pipeline recupera la señal sintética muy por encima del azar")
else:
    print("NOK! el pipeline no supera el azar")
