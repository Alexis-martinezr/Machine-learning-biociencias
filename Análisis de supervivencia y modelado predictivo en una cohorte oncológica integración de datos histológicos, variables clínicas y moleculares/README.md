**<h1>**🧬 Análisis de Supervivencia y Modelado Predictivo en una Cohorte Oncológica**</h1>**

Integración de datos histológicos, clínicos y moleculares

📘 **Introducción**

Este proyecto tiene como objetivo explorar, analizar y modelar la supervivencia en una cohorte oncológica mediante la integración de variables clínicas, histológicas y moleculares.
El dataset contiene información de 95 pacientes en una cohorte de validación, donde cada fila representa a un paciente individual.

📊 **Descripción del Dataset**
🧾 **Columnas principales**

- Patient ID: Identificador único (no se utiliza directamente en los modelos, pero es útil para trazabilidad).
- Survival time (days): Tiempo de seguimiento hasta el evento o censura (variable principal en análisis de supervivencia).
- Event: Indicador binario (1 = muerte/evento, 0 = vivo/censurado).
- Tumor size (cm): Tamaño del tumor (variable clínica continua).
- Grade (1–3): Diferenciación histológica (1 = bajo, 3 = alto, más agresivo).
- Stage (TNM 8th): Etapa clínica (IA, IB, II, III, IV, Recurrence).
- Age: Edad del paciente (años).
- Sex: Sexo biológico (Masculino / Femenino).
- Cigarette & Pack per year: Información sobre el hábito tabáquico.
- Type Adjuvant: Tipo de terapia adyuvante recibida (quimio, radio, etc.) – con valores ausentes.
- Batch: Lote de muestra o experimento (útil para corregir sesgos técnicos).
- EGFR / KRAS: Estado mutacional (Mutado / No mutado) – biomarcadores moleculares con valores faltantes.

🎯 **Objetivos del Proyecto**
El dataset permite realizar estudios en tres grandes líneas:
1. Análisis descriptivo y exploratorio
2. Distribuciones de edad, sexo, tabaquismo y mutaciones.
3. Comparación de supervivencia según estadio, grado, mutaciones o tipo de tratamiento.
4. Detección de patrones en pacientes con mayor o menor tiempo de supervivencia.
5. Análisis de supervivencia
6. Kaplan-Meier: Comparación de curvas de supervivencia entre grupos.
7. Regresión de Cox: Evaluación del efecto de variables (edad, estadio, mutaciones, tabaquismo) sobre la supervivencia.
8. Modelado predictivo (Machine Learning)
9. Predicción de supervivencia mediante modelos de clasificación supervisada.

🩺 **Clasificación TNM (8ª edición)**
- Estadio I (IA1, IA2, IA3, IB):
  Tumores pequeños sin ganglios ni metástasis. Generalmente curables con cirugía.

- Estadio II (IIA, IIB):
  Tumores más grandes o con invasión local. Tratamiento: cirugía + quimioterapia adyuvante.

- Estadio III (IIIA, IIIB):
  Tumores avanzados localmente o con afectación de ganglios mediastínicos. Tratamiento principal: radioquimioterapia.

- Estadio IV (IVA, IVB):
  Cáncer metastásico. Tratamiento: terapias sistémicas (quimio, inmunoterapia, terapias dirigidas).

- Recurrence:
  Reaparición del cáncer tras tratamiento inicial. No es un estadio, sino una categoría clínica.

🔬** Grado Histológico**
Describe la diferenciación celular del tumor:
- Grado 1: Bien diferenciado (crecimiento lento, mejor pronóstico).

- Grado 2: Moderadamente diferenciado (pronóstico intermedio).

- Grado 3: Poco diferenciado o indiferenciado (crecimiento rápido, peor pronóstico).

🤖 **Modelado Predictivo**

Se construyeron varios modelos de clasificación supervisada:

- Modelo	Accuracy (test)	Observaciones
- Árboles de decisión	~68%	Bajo poder predictivo
- KNN	~68%	Sensible a distribución y cantidad de datos
- Random Forest (5000 árboles)	AUC = 0.80	Mejores resultados

El modelo de Random Forest logró predecir correctamente el 85% de los verdaderos negativos y el 66.66% de los verdaderos positivos tras optimizar el threshold.

Las variables más relevantes fueron:
- Tumor size
- Grade
- Stage
- Sex
- Cigarette
- Pack per year
- EGFR
- KRAS
- Adjuvant.

⚠️ **Limitaciones:**

1. Tamaño pequeño del dataset.

2. Desbalance de clases (más negativos que positivos).

3. Posible falta de variables clave que expliquen la variabilidad.

📚 **Referencias**

- Mi, Haoyang (2023). Clinical_Data. Figshare. Dataset. https://doi.org/10.6084/m9.figshare.22638145.v1

🛠️ Tecnologías utilizadas

- Python (pandas, numpy, scikit-learn, lifelines, matplotlib, seaborn)

- Jupyter Notebook

- Git / GitHub

🚀 **Próximos pasos**

1. Recolectar más datos para mejorar el rendimiento del modelo.

2. Implementar técnicas de balanceo de clases (SMOTE, undersampling).

3. Probar modelos más avanzados (XGBoost, LightGBM, redes neuronales).

4. Validar el modelo en cohortes independientes.
