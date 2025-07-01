from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import os
import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem import rdFingerprintGenerator
import warnings
warnings.filterwarnings('ignore')

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Global variable to store the loaded model, scaler, and feature names
loaded_model_data = None

# Mapping from target protein names (from UI) to protein IDs (used in training)
# NOTE: 'AChE' is mapped to 'Q04844' to match the 5 one-hot encoded features from the training notebook.
PROTEIN_MAP = {
    'MAO-B': 'P27338',
    'COX-2': 'P35354',
    'VISFATIN': 'P43490',
    'BACE1': 'P56817',
    'AChE': 'Q04844'
}
# The order of columns must be consistent with the training data (get_dummies sorts them)
PROTEIN_IDS_ORDERED = sorted(PROTEIN_MAP.values())
TARGET_PROTEIN_NAMES = list(PROTEIN_MAP.keys())
DESCRIPTOR_NAMES = ["MolWt", "MolLogP", "NumHDonors", "NumHAcceptors", "TPSA"]

def load_model_data():
    """Load the serialized model dictionary (model, scaler, feature_names)"""
    global loaded_model_data
    try:
        pipeline_path = os.path.join(os.path.dirname(__file__), 'Model/final_tuned_model.pkl')
        if os.path.exists(pipeline_path):
            loaded_model_data = joblib.load(pipeline_path)
            print("Model data dictionary loaded successfully.")
            # Optional: Print keys to confirm content
            print(f"Loaded components: {list(loaded_model_data.keys())}")
            return True
        else:
            print(f"Model file 'final_tuned_model.pkl' not found at {pipeline_path}")
            return False
    except Exception as e:
        print(f"Error loading model data: {e}")
        return False

def smiles_to_features(smiles, target_protein_name):
    """
    Convert SMILES string and target protein to a feature vector that matches the training pipeline.
    Feature vector format: [Morgan Fingerprints (2048), One-Hot Proteins (5), Scaled Descriptors (5)]
    """
    try:
        # --- 1. Validate SMILES and Target Protein ---
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        if target_protein_name not in PROTEIN_MAP:
            raise ValueError(f"Invalid target protein name: {target_protein_name}")

        # --- 2. Generate Morgan Fingerprints (2048 bits) ---
        # CORRECTED: Using the more direct AllChem method
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
        fp_array = np.array(fp)

        # --- 3. One-Hot Encode Target Protein (5 features) ---
        target_protein_id = PROTEIN_MAP[target_protein_name]
        protein_array = np.zeros(len(PROTEIN_IDS_ORDERED))
        protein_array[PROTEIN_IDS_ORDERED.index(target_protein_id)] = 1

        # --- 4. Calculate Physicochemical Descriptors (5 features) ---
        descriptors = [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol)
        ]
        descriptor_df = pd.DataFrame([descriptors], columns=DESCRIPTOR_NAMES)
        
        # --- 5. Scale Descriptors using the loaded scaler ---
        if loaded_model_data and 'scaler' in loaded_model_data:
            scaled_descriptors = loaded_model_data['scaler'].transform(descriptor_df)
        else:
            raise RuntimeError("Feature scaler is not loaded.")
        
        scaled_descriptors_array = scaled_descriptors.flatten()

        # --- 6. Concatenate all features into a single vector (2058 features) ---
        final_features = np.concatenate([fp_array, protein_array, scaled_descriptors_array])
        
        # Reshape to (1, n_features) for the model
        return final_features.reshape(1, -1)

    except Exception as e:
        raise ValueError(f"Error processing SMILES '{smiles}': {str(e)}")


def make_prediction(smiles, target_protein):
    """
    Make a prediction for a given SMILES and target protein.
    """
    try:
        # Convert to features
        features = smiles_to_features(smiles, target_protein)

        # Make prediction using the loaded model
        if loaded_model_data and 'model' in loaded_model_data:
            model = loaded_model_data['model']
            # Ensure feature names match if the model requires it (e.g., some scikit-learn versions)
            # In this case, the model is a RandomForestClassifier which works on numpy arrays.
            prediction = model.predict(features)[0]
            prediction_proba = model.predict_proba(features)[0]
            confidence = max(prediction_proba)
        else:
            raise RuntimeError("Model is not loaded. Cannot make a prediction.")
        
        # Convert to human-readable format
        prediction_label = 'Active' if prediction == 1 else 'Inactive'
        
        return {
            'prediction': prediction_label,
            'confidence': float(confidence),
            'smiles': smiles,
            'target_protein': target_protein
        }
        
    except Exception as e:
        raise ValueError(f"Prediction error: {str(e)}")


@app.route('/')
def index():
    """Serve the main page"""
    return render_template('index.html')


@app.route('/predict', methods=['POST'])
def predict():
    """
    Handle single compound prediction. Supports both single target and 'All Targets' mode.
    """
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        smiles = data.get('smiles', '').strip()
        target_protein = data.get('target_protein', '').strip()
        
        if not smiles:
            return jsonify({'error': 'SMILES string is required'}), 400
        if not target_protein:
            return jsonify({'error': 'Target protein is required'}), 400
        
        # Handle 'All Targets' mode
        if target_protein == 'All Targets':
            results = []
            for protein_name in TARGET_PROTEIN_NAMES:
                try:
                    result = make_prediction(smiles, protein_name)
                    results.append({
                        'target': protein_name,
                        'prediction': result['prediction'],
                        'confidence': result['confidence']
                    })
                except Exception as e:
                    results.append({
                        'target': protein_name,
                        'prediction': 'Error',
                        'confidence': 0.0,
                        'error': str(e)
                    })
            return jsonify(results)
        
        # Handle single target prediction
        else:
            result = make_prediction(smiles, target_protein)
            return jsonify(result)
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/predict_batch', methods=['POST'])
def predict_batch():
    """
    Handle batch prediction. Expects a list of {smiles: "...", target_protein: "..."} objects.
    """
    try:
        data = request.get_json()
        if not data or not isinstance(data, list):
            return jsonify({'error': 'Expected a list of prediction requests'}), 400
        
        results = []
        for item in data:
            if not isinstance(item, dict):
                results.append({'error': 'Each item must be a dictionary with smiles and target_protein'})
                continue
            
            smiles = item.get('smiles', '').strip()
            target_protein = item.get('target_protein', '').strip()
            
            if not smiles or not target_protein:
                results.append({
                    'smiles': smiles, 'target_protein': target_protein,
                    'error': 'Both smiles and target_protein are required'
                })
                continue
            
            try:
                result = make_prediction(smiles, target_protein)
                results.append(result)
            except Exception as e:
                results.append({
                    'smiles': smiles, 'target_protein': target_protein,
                    'error': str(e)
                })
        
        return jsonify(results)
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    model_loaded = loaded_model_data is not None
    return jsonify({
        'status': 'healthy' if model_loaded else 'unhealthy',
        'model_loaded': model_loaded,
        'available_targets': TARGET_PROTEIN_NAMES
    })


if __name__ == '__main__':
    # Load the model data on startup
    load_model_data()
    
    # Run the app
    app.run(host='0.0.0.0', port=5000, debug=True)
