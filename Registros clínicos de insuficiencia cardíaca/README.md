**<h1>**Análisis de Registros Clínicos de Insuficiencia Cardíaca ❤️ **</h1>**

Este proyecto analiza datos clínicos de pacientes con insuficiencia cardíaca, con el objetivo de construir modelos predictivos de supervivencia basados en variables clínicas clave.

📘 **Introducción**

El análisis se fundamenta en los datos publicados por Chicco y Jurman (2020), correspondientes a una cohorte de 299 pacientes registrados entre abril y diciembre de 2015 en hospitales de Faisalabad, Pakistán (BioMed Central).

El objetivo original del estudio fue comparar distintos clasificadores para predecir la supervivencia, encontrando que solo dos variables clínicas —creatinina sérica y fracción de eyección— eran suficientes para construir modelos predictivos precisos.

En este proyecto se desarrollaron nueve modelos basados en:

- Bosques Aleatorios (Random Forest)

- Regresión Logística

- SVC (Support Vector Classifier)

📊 **Distribución Demográfica**

- Sexo: 105 mujeres, 194 hombres

- Edad: 40 a 95 años

Condición clínica: **Todos con disfunción sistólica del ventrículo izquierdo, clases III o IV según NYHA**

🧾 **Variables Incluidas**

Categóricas (binarias): 
- Anemia
- Hipertensión
- Diabetes
- Sexo
- Tabaquismo

Continuas: 
- Creatinina fosfoquinasa (CPK)
- Fracción de eyección
- Creatinina sérica
- Sodio sérico

Variable objetivo: **death_event (1 = falleció, 0 = sobrevivió)**

⚖️ **Balance de clases**

96 fallecimientos (~32%)

203 sobrevivientes (~68%)

🎯 **Objetivos del Proyecto**

1. Validar los hallazgos de Chicco y Jurman respecto a las variables más determinantes para predecir supervivencia.

2. Construir modelos predictivos de supervivencia usando diferentes técnicas de machine learning.

3. Optimizar el threshold para balancear sensibilidad y especificidad según necesidades clínicas.

3. Evaluar el impacto de variables adicionales (edad, sodio sérico, hipertensión, anemia, CPK, plaquetas, tabaquismo, sexo, diabetes) en la predicción.

🤖 **Modelado Predictivo**

El mejor desempeño se obtuvo con un modelo de Random Forest, alcanzando:

- AUC ROC: 0.80 (sin tiempo de seguimiento)

Tras ajustar el **threshold a 0.597**, se logró una **precisión del 81.1%**, balanceando verdaderos positivos y negativos.

Al incluir el tiempo de seguimiento, el **AUC ROC** aumentó hasta **0.90**.

**Las variables más influyentes:** creatinina sérica, fracción de eyección y edad.

Variables adicionales como sodio sérico, hipertensión, anemia, CPK, plaquetas, tabaquismo, sexo y diabetes contribuyeron a mejorar el rendimiento en combinación con el tiempo de seguimiento.

🛠️ **Tecnologías Utilizadas**

- Python: pandas, numpy, scikit-learn, matplotlib, seaborn

- Jupyter Notebook

- Git / GitHub

📚 **Referencias**

- Chicco, D., Jurman, G. Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. BMC Med Inform Decis Mak 20, 16 (2020). https://doi.org/10.1186/s12911-020-1023-5

🚀 **Próximos pasos**

1. Explorar técnicas de balanceo de clases (SMOTE, undersampling) para mejorar detección de fallecimientos.

2. Probar modelos avanzados como XGBoost o LightGBM.

3. Validar el modelo en cohortes independientes.

4. Analizar la relevancia clínica de otras variables adicionales para mejorar la predicción.
