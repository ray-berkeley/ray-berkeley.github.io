#!/usr/bin/env python3
import cgi
import cgitb
import os
import sys
import json
import tempfile
import traceback

# Enable CGI traceback for debugging
cgitb.enable()

# Import your MyCCD module (ensure myccd.py is in the same directory)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from myccd import MyCCD, generate_ccd_from_smiles, generate_ccd_from_pdb

def main():
    # Print CGI header
    print("Content-Type: application/json")
    print()  # Empty line to separate headers from body

    try:
        # Parse form data
        form = cgi.FieldStorage()

        # Get form fields
        name = form.getvalue('name', 'PROTEIN-LIGAND')
        protein_sequence = form.getvalue('protein_sequence', '')
        dna_sequences_text = form.getvalue('dna_sequences', '')
        ligand_input_type = form.getvalue('ligand_input_type', 'smiles')
        ligand_smiles = form.getvalue('ligand_smiles', '')
        ligand_name = form.getvalue('ligand_name', 'LIG')
        protein_residue = int(form.getvalue('protein_residue', 0))
        protein_atom = form.getvalue('protein_atom', 'SG')
        ligand_atom = form.getvalue('ligand_atom', 'C1')
        model_seeds_text = form.getvalue('model_seeds', '42')

        # Parse model seeds
        model_seeds = [int(seed.strip()) for seed in model_seeds_text.split(',')]

        # Process DNA sequences
        dna_sequences = []
        for i, seq in enumerate(dna_sequences_text.strip().split('\n')):
            if seq.strip():
                dna_sequences.append({
                    "dna": {
                        "id": chr(66 + i),  # B, C, D...
                        "sequence": seq.strip()
                    }
                })

        # Process ligand
        ccd = None

        if ligand_input_type == 'smiles' and ligand_smiles:
            # Process SMILES string
            ccd = generate_ccd_from_smiles(ligand_smiles, ligand_name)
        else:
            # Process uploaded file
            if 'ligand_file' in form:
                ligand_file = form['ligand_file']
                if hasattr(ligand_file, 'filename') and ligand_file.filename:
                    # Save uploaded file temporarily
                    temp_fd, temp_path = tempfile.mkstemp(suffix='_' + os.path.basename(ligand_file.filename))
                    os.close(temp_fd)

                    with open(temp_path, 'wb') as f:
                        f.write(ligand_file.file.read())

                    # Determine file type by extension
                    file_ext = os.path.splitext(ligand_file.filename)[1].lower()

                    if file_ext == '.pdb':
                        ccd = generate_ccd_from_pdb(temp_path, ligand_name)
                    else:
                        # Assume SMILES file or direct SMILES string
                        ccd = generate_ccd_from_smiles(temp_path, ligand_name)

                    # Clean up temp file
                    os.unlink(temp_path)

        if not ccd:
            error_result = {"error": "No valid ligand input provided"}
            print(json.dumps(error_result))
            return

        # Generate CCD content (mmCIF format)
        temp_cif_fd, temp_cif_path = tempfile.mkstemp(suffix='.cif')
        os.close(temp_cif_fd)
        ccd.write_ccd(temp_cif_path)

        with open(temp_cif_path, 'r') as f:
            ccd_content = f.read()

        # Clean up temp CIF file
        os.unlink(temp_cif_path)

        # Prepare JSON structure
        json_output = {
            "name": name,
            "modelSeeds": model_seeds,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": protein_sequence.strip()
                    }
                }
            ],
            "bondedAtomPairs": [
                [
                    ["A", protein_residue, protein_atom],
                    ["D", 1, ligand_atom]
                ]
            ],
            "userCCD": ccd_content,
            "dialect": "alphafold3",
            "version": 2
        }

        # Add DNA sequences
        json_output["sequences"].extend(dna_sequences)

        # Add ligand entry
        json_output["sequences"].append({
            "ligand": {
                "id": "D",
                "ccdCodes": [ligand_name]
            }
        })

        # Output JSON response
        print(json.dumps(json_output, indent=2))

    except Exception as e:
        # Log the full traceback
        error_traceback = traceback.format_exc()
        sys.stderr.write(error_traceback)

        # Return error JSON
        error_result = {"error": str(e)}
        print(json.dumps(error_result))

if __name__ == '__main__':
    main()