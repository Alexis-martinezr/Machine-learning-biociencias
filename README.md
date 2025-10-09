**<h1>**📂 **Portafolio de Proyectos de Machine Learning y Análisis Clínico****</h1>**

Este repositorio contiene cuatro proyectos de análisis de datos clínicos y modelado predictivo utilizando técnicas de Machine Learning, Bioinformática y análisis estadístico. Cada proyecto está organizado en su propia carpeta, con un README.md detallado por proyecto.

📝 **Proyectos incluidos**
1️⃣ **Análisis de Supervivencia en Cohorte Oncológica**

- Carpeta: cancer_survival_analysis

- Descripción: Estudio de supervivencia en 95 pacientes oncológicos, integrando datos histológicos, clínicos y moleculares.

- Técnicas: Kaplan-Meier, regresión de Cox, modelos de árboles de decisión, KNN y Random Forest.

- Resultados clave: Mejor modelo de Random Forest con AUC 0.8 y optimización de threshold, destacando la importancia de variables como tamaño del tumor, grado, estadio, EGFR y KRAS.

2️⃣ **Análisis de Registros Clínicos de Insuficiencia Cardíaca**

- Carpeta: heart_failure_analysis

- Descripción: Análisis de 299 pacientes con insuficiencia cardíaca para predecir eventos de mortalidad usando variables clínicas.

- Técnicas: Random Forest, Regresión Logística, SVC, optimización de thresholds.

- Resultados clave: Creatinina sérica y fracción de eyección son las variables más determinantes; mejor modelo Random Forest AUC 0.8–0.9 según inclusión de tiempo de seguimiento.

3️⃣ **Optimización de Modelos de Detección Temprana de Cáncer de Mama**

- Carpeta: breast_cancer_knn

- Descripción: Clasificación de tumores benignos vs. malignos usando el dataset load_breast_cancer de scikit-learn.

- Técnicas: K-Nearest Neighbors (KNN), PCA, optimización de hiperparámetros, análisis de threshold.

- Resultados clave: La reducción de dimensionalidad y normalización mejoran el desempeño del modelo; selección de threshold entre 0.1–0.3 maximiza la detección de verdaderos positivos.

4️⃣ **Análisis de Datos de COVID-19 en México (IMSS) (en progreso)**

- Carpeta: covid_imss_analysis

- Descripción: Análisis exploratorio de datos de pacientes COVID-19 en México, provenientes del IMSS.

- Técnicas: Limpieza y visualización de datos, análisis descriptivo de variables clínicas y demográficas.

- Estado: Proyecto en desarrollo; se planea incluir modelado predictivo y comparaciones temporales.

🛠️ **Tecnologías utilizadas**

- Python (pandas, numpy, scikit-learn, matplotlib, seaborn)

- Jupyter Notebook

- Git / GitHub

📚 **Referencias**

- Chicco, D., Jurman, G. Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. BMC Med Inform Decis Mak 20, 16 (2020). https://doi.org/10.1186/s12911-020-1023-5

- Mi, Haoyang (2023). Clinical_Data. figshare. https://doi.org/10.6084/m9.figshare.22638145.v1

- scikit-learn: load_breast_cancer dataset

- IMSS COVID-19 Dataset (México)

🚀 **Próximos pasos**

1. Validar modelos en cohortes externas o datasets independientes.

2. Probar modelos avanzados (XGBoost, LightGBM) y técnicas de balanceo de clases.

3. Mejorar la documentación y visualizaciones de los proyectos.

4. Completar y expandir el análisis de COVID-19 con modelado predictivo y comparaciones temporales.
