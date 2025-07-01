<p align="center">
  <img alt="UI" src="https://i.imgur.com/irO3zIU.png">
</p>

### 🌐 Overview

Welcome to our cutting-edge computational pipeline designed to accelerate Alzheimer's Disease (AD) research. This project integrates advanced bioinformatics and cheminformatics, creating a seamless workflow from raw **single-cell RNA sequencing (scRNA-seq)** data to predictive **Quantitative Structure-Activity Relationship (QSAR)** modeling.

Our mission is to democratize access to powerful predictive tools, lowering the barrier to entry for researchers in the neurodegenerative disease space. This repository provides a comprehensive toolkit for data integration, cellular analysis, and machine learning-based bioactivity prediction.

You can access and use the live application at: https://QSARify.com

<!-- Placeholder for a high-level architecture diagram -->
<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

### 🔗 Project Workflow
<p align="center">
  <img alt="UI" src="https://i.imgur.com/9aMZA6x.png">
</p>

<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

### ✨ Features

This pipeline is organized into three core modules, each providing a distinct set of functionalities.

### 🔬 Single-Cell Analysis
* **🔗 Data Integration**: Merges complex `scRNA-seq` datasets from multiple public studies (**GSE138852**, **GSE157827**, **GSE175814**, **GSE163577**) across various brain regions into a unified `Seurat` object.
* * Notably, the original studies encompassed over 500,000 cells, all of which were processed in our analysis. For ease of GitHub upload, a random subset of 25,000 cells (approximately 5,000 from each study) has been provided
* **✅ Quality Control & Normalization**: Implements a robust QC pipeline that filters low-quality cells (e.g., mitochondrial content > 10%), normalizes data using `SCTransform`, and removes artifacts with `DoubletFinder`.
* **🧬 Cell Type Annotation**: Automatically annotates cell clusters using the `SingleR` package with established reference datasets.
* **📊 Dimensional Reduction & Clustering**: Leverages `PCA` for initial dimensionality reduction, `Harmony` for batch effect correction, and `UMAP` for visualization and clustering.
* **🔍 Differential Abundance**: Employs `MiloR` to identify statistically significant changes in cell population abundance between experimental conditions.
* **📡 Cell-Cell Communication**: Uses `CellChat` to infer and analyze intercellular communication networks, identifying key ligand-receptor interactions and signaling pathways.

### 🧪 QSAR Modeling & Bioactivity Prediction
* **🎯 Target Data Collection**: Queries the **ChEMBL** database using UniProt IDs for key AD-related targets (`MAO-B`, `COX-2`, `VISFATIN`, `BACE1`, `AChE`). It retrieves bioactivity data (e.g., *IC50* values) and calculates critical **ADME** properties (*MW*, *LogP*, *HBD*, *HBA*, *Lipinski's Rule*).
* **🧠 Machine Learning Pipeline**:
    * **Feature Engineering**: Converts **SMILES** strings into 2048-bit **Morgan fingerprints** and calculates key physicochemical descriptors.
    * **Imbalance Handling**: Utilizes the `SMOTETomek` technique to address class imbalance in the bioactivity data.
    * **Model Training & Tuning**: Trains multiple baseline models (*Random Forest*, *Gradient Boosting*, *Neural Network*) and hyperparameter-tunes the top performer (*Random Forest*) to achieve an **AUC of ~0.82** on the test set.
* **📈 Model Evaluation**: Generates comprehensive visualizations for model performance, including ROC curves, accuracy plots, and F1-score comparisons.
    <!-- Placeholder for the ROC curve image -->
    <p align="center">
    </p>

### 🌐 Interactive Web API
* **🚀 RESTful Endpoints**: A `Flask`-based API provides endpoints for health checks (`/health`), single predictions (`/predict`), and batch predictions (`/predict_batch`).
* **🎨 User-Friendly Interface**: An intuitive web UI for interacting with the QSAR model.
    * **Single & Batch Predictions**: Supports bioactivity prediction for one or multiple compounds via SMILES input or file upload (`.txt`, `.xls`, `.xlsx`).
    * **Dynamic Target Selection**: Allows users to predict against a specific target or all available targets.
    * **History Tracking**: A sortable and filterable history tab keeps a record of all prediction tasks.
    * **About Page**: Contains project details and team information.
 
<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

<p align="center">
  <img alt="UI" src="https://i.imgur.com/C0RI7QT.png">
</p>

<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

## 📂 Project Structure
```
AlzheimerDisease_FromSingleCell/
├── SingleCell/                # 🧬 scRNA-seq preprocessing (R)
│   ├── Merge_Data.R
│   └── SingleCell_Main.R
├── MiloR/                     # 🧬 Differential-abundance analysis (R)
│   └── MiloR_CellAbundance.R
├── CellChat/                  # 🧬 Cell–cell communication analysis (R)
│   └── CellChat.R
├── QSAR/                      # 🧠 QSAR modeling & web app (Python)
│   ├── figures/
│   │   └── roc_curves_comparison.png
│   ├── Data/
│   │   ├── chembl_results_P_27338_MAO-B_IC50_classified.csv
│   │   ├── chembl_results_P_35354_COX2_IC50_classified.csv
│   │   ├── chembl_results_P_43490_VISFATIN_IC50_classified.csv
│   │   ├── chembl_results_P_56817_BACE1_IC50_classified.csv
│   │   └── chembl_results_Q_04844_ACHE_IC50_classified.csv
│   ├── Model/
│   │   └── final_tuned_model.pkl
│   ├── templates/
│   │   └── index.html
│   ├── Target_Collection.ipynb
│   ├── Ligand_Final.ipynb
│   ├── app.py
│   └── requirements.txt
├── Data/                      # A large-scale analysis of over 500,000 cells was performed. A 25,000-cell subset (5,000 from each study) is provided on GitHub for convenience.
│   └── 25K_Sample.rds
└── README.md
```

<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

## 📝 Detailed Script Information

| Script 🖥️                        | Purpose 🎯                                                                   | Key Libraries 🛠️                                                              | Output 📄                                                            |
| -------------------------------- | ---------------------------------------------------------------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| `Merge_Data.R`                   | Integrates raw scRNA-seq count matrices from multiple GSE studies.             | `Seurat`, `batchelor`, `SingleCellExperiment`                                 | A unified Seurat object containing all datasets.                     |
| `SingleCell_Main.R`              | Performs QC, normalization, clustering, and cell type annotation.            | `Seurat`, `harmony`, `DoubletFinder`, `SingleR`                               | A processed Seurat object with UMAPs and cell annotations.           |
| `MiloR_CellAbundance.R`       | Conducts differential abundance testing on cell neighborhoods.               | `miloR`, `SingleCellExperiment`, `ggplot2`                                    | Differential abundance statistics and visualizations.                |
| `CellChat.R`  | Infers and analyzes cell-cell communication pathways.                        | `CellChat`, `Seurat`, `dplyr`                                                 | Communication network data and plots (bubble plots, heatmaps).       |
| `Target_Collection.ipynb`        | Retrieves and preprocesses bioactivity & ADME data from ChEMBL.              | `pandas`, `chembl_webresource_client`, `rdkit`                                | A cleaned DataFrame and exploratory data visualizations.             |
| `Ligand_Final.ipynb`             | Trains, tunes, and evaluates the QSAR machine learning model.                | `scikit-learn`, `imbalanced-learn`, `rdkit`, `pandas`                         | A serialized model (`.pkl`) and performance plots.                   |
| `app.py`                         | Serves a Flask-based web API for on-demand bioactivity predictions.          | `flask`, `flask-cors`, `joblib`, `rdkit`                                      | JSON responses with predictions and confidence scores.               |
| `index.html`                     | Provides the interactive front-end UI for the QSAR prediction tool.          | HTML, CSS, JavaScript                                                         | An interactive web interface rendered in the browser.                |

<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

## 🛠️ Installation

To set up the project environment, please follow these steps.

#### **Step 1: Install R Dependencies** 📦
(Required for `SingleCell`, `MiloR`, and `CellChat` analysis)
```R
# Install core packages from CRAN
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "scater", "scran", "harmony", "batchelor", "SingleR"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "miloR", "glmGamPoi"))
```

#### **Step 2: Install Python Dependencies** 🐍
(Required for `QSAR` modeling and the `Flask` API)
```bash
# Clone the repository
git clone [https://github.com/xhammady/AD-scRNA2QSAR.git](https://github.com/xhammady/AD-scRNA2QSAR.git)
cd AD-scRNA2QSAR

# Install Python packages from requirements.txt
pip install -r QSAR/requirements.txt
```
> *Note: Key Python libraries include `chembl-webresource-client`, `rdkit`, `scikit-learn`, `imbalanced-learn`, `pandas`, `flask`, and `flask-cors`.*

<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

## 🚀 Usage

Follow this sequence to run the full analysis pipeline.

#### **➡️ Step 1: Run the Single-Cell Bioinformatics Pipeline**
1.  **🧬 Merge Data**: Execute `SingleCell/Merge_Data.R` to combine the raw count matrices.
2.  **🔬 Preprocess & Cluster**: Run `SingleCell/SingleCell_Main.R` to perform QC, normalization, integration, and annotation.
3.  **🔍 Analyze Differential Abundance**: Use `MiloR/MiloR_CellAbundance.R` to compare cell populations.
4.  **📡 Infer Communication**: Run `CellChat/CellChat.R` to analyze signaling pathways.

#### **➡️ Step 2: Run the QSAR Cheminformatics Pipeline**
1.  **🧪 Collect Target Data**: Open and run the `QSAR/Target_Collection.ipynb` notebook to query ChEMBL and generate the analysis dataset.
2.  **🧠 Train ML Model**: Open and run `QSAR/Ligand_Final.ipynb` to preprocess features, train the Random Forest model, and save the final `.pkl` file.

#### **➡️ Step 3: Launch the Predictive API & User Interface**
1.  **🌐 Start the Server**: From the command line, run the Flask application:
    ```bash
    python QSAR/app.py
    ```
2.  🎨 Access the UI: Open your web browser and navigate to http://localhost:5000/ or visit the live application at https://QSARify.com. You can now:
    * Enter a **SMILES** string for a single compound prediction.
    * Upload a file (`.txt`, `.xls`, `.xlsx`) for batch predictions.
    * View and manage results in the **History** tab.

<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>

## 🤝 Contributing
We welcome contributions to improve this project! Please fork the repository, create a new branch for your feature, and submit a pull request with a detailed description of your changes. Ensure you follow existing coding standards and include tests where applicable.

## 📜 License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## 🙏 Acknowledgments & Contributors

#### **Project Team**
* **Amr Mohamed ElHefnawy**
* **Ibrahim Abdelkarim Hammad**
* **Mira Moheb Attia**
* **Abdelrahman Wagih**
* **Reem Sharaf EL-Deen Hassan**
* **Lorance Gergis Labeeb**
* **Prof. Dr. Sameh Ibrahim Hassanien**

#### **Tools & Data**
* **🧬 Software & Libraries**: `Seurat`, `CellChat`, `MiloR`, `RDKit`, `scikit-learn`, `imbalanced-learn`, `Flask`, `DoubletFinder`, `SingleR`, `pandas`, `numpy`, `matplotlib`.
* **🧪 Databases**: We gratefully acknowledge the **ChEMBL** database for providing essential bioactivity data.
* **📊 Data Providers**: This work would not be possible without the public datasets provided by Gene Expression Omnibus (GEO): **GSE138852**, **GSE157827**, **GSE175814**, **GSE163577**.

## 📧 Contact
* For questions or issues, please open an issue on GitHub or contact the project maintainers through mail.
* [**Amr Mohamed ElHefnawy**](mailto:AAlhfnawy@nu.edu.eg)
* [**Ibrahim Abdelkarim Hammad**](mailto:i.abdelkarim@nu.edu.eg)
* [**Mira Moheb Attia**](mailto:m.moheb2119@nu.edu.eg)
* [**Abdelrahman Wagih**](mailto:a.wagih2160@nu.edu.eg)
* [**Reem Sharaf EL-Deen Hassan**](mailto:r.sharaf2187@nu.edu.eg)
* [**Lorance Gergis Labeeb**](mailto:l.gergis2184@nu.edu.eg)
* [**Prof. Dr. Sameh Ibrahim Hassanien**](mailto:SIbrahem@nu.edu.eg)
<p align="center">
  <img alt="border" src="https://i.imgur.com/obXZbsy.png">
</p>
