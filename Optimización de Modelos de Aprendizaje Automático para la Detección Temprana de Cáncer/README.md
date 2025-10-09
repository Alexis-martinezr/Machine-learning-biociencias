**<h1>**🩺 **Optimización de Modelos de Machine Learning para la Detección Temprana de Cáncer de Mama****</h1>**

Este proyecto formativo se centra en la construcción y optimización de modelos de machine learning para clasificar tumores como malignos o benignos utilizando el dataset load_breast_cancer de scikit-learn.

📘 **Introducción**

El análisis se realizó con el objetivo de:

Evaluar la correlación entre las variables independientes y la variable dependiente (diagnóstico de tumor).

Construir modelos K-Nearest Neighbors (KNN) para la clasificación de tumores.

Analizar el efecto de la reducción de dimensionalidad y la normalización de datos sobre el desempeño del modelo.

Debido a que varias variables independientes presentaron alta correlación, se implementó reducción de dimensionalidad mediante PCA, seleccionando la cantidad óptima de componentes para maximizar la precisión del modelo.

📊 **Dataset**

Fuente: **load_breast_cancer de scikit-learn**

Número de muestras: **569 pacientes**

Variables independientes: **30 características clínicas y bioquímicas**

**Variable objetivo:** target

- **0 = Benigno**

- **1 = Maligno**

**Observaciones importantes**

La variable objetivo presenta correlación negativa con la mayoría de las variables independientes, lo que indica que ciertos valores altos de estas características se asocian con un tumor maligno.

Algunas variables independientes muestran alta correlación entre sí, lo que motivó la construcción de modelos con y sin reducción de dimensionalidad.

🎯 **Objetivos del Proyecto**

1. Construir modelos KNN para la clasificación de tumores.

2. Comparar modelos con y sin reducción de dimensionalidad (PCA) y normalización de parámetros.

3. Optimizar el umbral de decisión (threshold) para maximizar la detección de verdaderos positivos (sensibilidad), minimizando falsos negativos.

4. Analizar el impacto de la selección de componentes y la normalización en el desempeño del modelo.

🤖 **Modelado y Optimización**

Se construyeron dos modelos principales:

1. Con reducción de dimensionalidad y normalización

2. Sin reducción de dimensionalidad ni normalización

La selección del threshold óptimo se realizó evaluando valores entre 0.1 y 0.3:

Thresholds de 0.2 y 0.3 lograron 8 falsos positivos, mientras que 0.1 produjo 10 falsos positivos.

La selección final del umbral debería ser realizada por un profesional médico, priorizando la detección de verdaderos positivos incluso si se incrementan los falsos positivos.

🛠️ *Tecnologías utilizadas*

- Python (scikit-learn, pandas, numpy, matplotlib, seaborn)

- Jupyter Notebook

- Git / GitHub

📚 **Referencias**

- Chicco, D. & Jurman, G. Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. BMC Med Inform Decis Mak 20, 16 (2020).

scikit-learn: load_breast_cancer dataset

🚀 **Próximos pasos**

1. Evaluar otros clasificadores (Random Forest, SVM, Logistic Regression) y compararlos con KNN.

2. Implementar técnicas de balanceo de clases si fuese necesario.

3. Optimizar aún más los hiperparámetros de KNN y PCA.

4. Validar modelos en datasets externos para confirmar la generalización.
