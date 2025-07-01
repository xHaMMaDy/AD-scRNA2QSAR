# Alzheimer's Drug Screening Tool 🧬

A professional machine learning-powered web application designed to predict the bioactivity of chemical compounds against key Alzheimer's disease-related protein targets.

## Project Overview 🌟

This project delivers a robust computational tool to support researchers in evaluating chemical compounds for bioactivity against critical Alzheimer's disease protein targets: MAO-B, COX-2, VISFATIN, BACE1, and AChE. Built with a RandomForestClassifier trained on advanced cheminformatics features—including Morgan fingerprints, physicochemical descriptors, and one-hot encoded protein targets—this application enables rapid predictions of a compound’s "Active" or "Inactive" status using SMILES string inputs. It offers both single compound and batch prediction capabilities through an intuitive, user-friendly interface, accelerating early-stage drug discovery.

## Features ✨

- **Single Compound Prediction**: 📥 Enter a SMILES string and select a target protein (or "All Targets" for a comprehensive disease panel) to receive bioactivity predictions with confidence scores.
- **Batch Prediction**: 📊 Process multiple SMILES strings against a single protein target, with support for text input or file uploads (.txt, .xls, .xlsx).
- **Interactive Web Interface**: 🌐 Developed using HTML, JavaScript, and Flask for a seamless and accessible user experience.
- **Advanced Machine Learning Model**: ⚙️ Powered by a pre-trained RandomForestClassifier leveraging Morgan fingerprints (2048 bits), one-hot encoded protein features, and scaled physicochemical descriptors (MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA).
- **Health Check Endpoint**: 🛠️ API endpoint (`/health`) to confirm model availability and list supported protein targets.
- **Scalable Backend**: 🔧 Flask-based server with CORS support for efficient front-end integration.

## Tech Stack 🛠️

- **Backend**: Python, Flask, RDKit, scikit-learn, pandas, NumPy, joblib
- **Frontend**: HTML, JavaScript, CSS
- **Machine Learning**: RandomForestClassifier, StandardScaler
- **Data Processing**: Morgan fingerprints (RDKit), physicochemical descriptors, one-hot encoding
- **Development Environment**: Jupyter Notebook for model training and feature engineering

## Installation ⚙️

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/xhammady/alzheimers-drug-screening.git
   cd alzheimers-drug-screening
   ```

2. **Set Up a Virtual Environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Download Pre-trained Model**:
   - Place the `final_tuned_model.pkl` file (containing the model, scaler, and feature names) in the project root directory.
   - Note: Due to size constraints, the model file is not included in the repository. Contact the team or use the training notebook (`Ligand_Final_Latest_Feature.ipynb`) to generate it.

5. **Run the Flask Application**:
   ```bash
   python app.py
   ```
   - Access the application at `http://localhost:5001`.

## Usage 📖

1. **Access the Web Interface**:
   - Navigate to `http://localhost:5001` in your browser.
   - Use the **Single Compound** tab to input a SMILES string (e.g., `C#CCN(C)[C@H](C)Cc1ccccc1`) and select a target protein.
   - Use the **Batch Prediction** tab to input multiple SMILES strings or upload a supported file format.

2. **API Endpoints**:
   - **Single Prediction**: `POST /predict`
     ```json
     {
       "smiles": "C#CCN(C)[C@H](C)Cc1ccccc1",
       "target_protein": "MAO-B"
     }
     ```
   - **Batch Prediction**: `POST /predict_batch`
     ```json
     [
       {"smiles": "C#CCN(C)[C@H](C)Cc1ccccc1", "target_protein": "MAO-B"},
       {"smiles": "O=C(O)c1coc2ccccc2c1=O", "target_protein": "MAO-B"}
     ]
     ```
   - **Health Check**: `GET /health`
     Returns model status and available protein targets.

3. **Training the Model**:
   - Execute the `Ligand_Final_Latest_Feature.ipynb` notebook to preprocess data, train the RandomForestClassifier, and save the model as `final_tuned_model.pkl`.
   - Ensure ChemBL dataset CSV files (`chembl_results_*.csv`) are available in the notebook’s directory.

## Project Structure 📂

```plaintext
├── app.py                    # Flask application backend
├── index.html                # Frontend web interface
├── Ligand_Final_Latest_Feature.ipynb  # Jupyter notebook for model training
├── final_tuned_model.pkl     # Pre-trained model (not included)
├── requirements.txt          # Python dependencies
└── README.md                   # Project documentation
```

## Model Details 📈

- **Input Features** (2058 total):
  - Morgan fingerprints: 2048 bits (radius=3)
  - One-hot encoded protein IDs: 5 features (P27338, P35354, P43490, P56817, Q04844)
  - Scaled physicochemical descriptors: 5 features (MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA)
- **Model**: RandomForestClassifier (n_estimators=100, random_state=42)
- **Preprocessing**: StandardScaler for descriptors, one-hot encoding for proteins
- **Training Data**: Combined ChemBL datasets for five protein targets (17,425 samples)
- **Performance**: Achieves ~86.83% accuracy and 0.8203 AUC-ROC on the test set.

## Contributing 🤝

We welcome contributions from the community! To contribute:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/your-feature`).
3. Commit your changes (`git commit -m "Add your message"`).
4. Push to the branch (`git push origin feature/your-feature`).
5. Open a pull request.

Please adhere to PEP 8 coding standards and include relevant tests with your submissions.

## Team 👥

- **Ibrahim Abdelkarim Hammad** (202001284)
- **Mira Moheb Attia** (211002019)
- **Abdelrahman Wagih** (211002381)
- **Reem Sharaf EL-Deen Hassan** (211001887)
- **Lorance Gergis Labeeb** (211002084)
- **Teaching Assistant**: Mr. Amr Mohamed ElHefnawy
- **Professor**: Dr. Sameh Ibrahim Hassanien

## License 📜

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments 🙏

- **ChemBL**: For providing high-quality bioactivity datasets.
- **RDKit**: For powerful cheminformatics tools.
- **scikit-learn**: For robust machine learning utilities.
- **Flask**: For enabling efficient web application development.