# TFM pipeline de integración de datos ómicos de cáncer de mama
# puesta en marcha

```bash
cd src
docker compose up -d
python manage.py migrate
python manage.py load_hgnc_local
python manage.py populate_pam50
python manage.py reimport_local
python manage.py sync_geo
```

# analisis de clasificacion

```bash
jupyter notebook notebooks/ml_pam50_classifier.ipynb # clasific interna TCGA-BRCA
jupyter notebook notebooks/cross_validation.ipynb # validac cross-platform GEO GSE25055
```

# memoria (LaTeX)

La memoria se compila con XeLaTeX desde `tex/`:

```bash
cd tex
make
```

La plantilla usa la fuente **Calibri**, que es propietaria de Microsoft y por eso no se redistribuye en este repositorio. Para compilar, copia los ficheros de Calibri incluidos en Windows/Office (`calibri.ttf`, `calibrib.ttf`, `calibrii.ttf`, `calibriz.ttf`, `calibril.ttf` y `calibrili.ttf`) a la carpeta `tex/fonts/`.