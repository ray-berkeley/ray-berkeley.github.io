#!/usr/bin/env python3
import os
import sys
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt

class MyCCD:
    """Class to represent a CDD-formatted mmCIF file."""

    def __init__(self, name=None, mol=None):
        """Initialize the MyCCD object.

        Args:
            name (str, optional): The name of the compound. Defaults to None.
            mol (rdkit.Chem.rdchem.Mol, optional): RDKit molecule object. Defaults to None.
        """
        self.name = name or "UNK"
        self.iupac_name = None
        self.formula = None
        self.weight = None
        self.pdb = None
        self.smiles = None
        self.mol = mol

        # Generate formula and weight from molecule if available
        if self.mol is not None:
            self.formula = CalcMolFormula(self.mol)
            self.weight = CalcExactMolWt(self.mol)

            # Generate CCD blocks
            self.comp = self.mol_to_ccd_comp_block(self.mol, self.name)
            self.comp_atom = self.mol_to_ccd_comp_atom_block(self.mol, self.name)
            self.comp_bond = self.mol_to_ccd_comp_bond_block(self.mol, self.name)
        else:
            self.comp = None
            self.comp_atom = None
            self.comp_bond = None

    def __str__(self):
        """Return a string representation of the MyCCD object."""
        return f"MyCCD object: {self.name}"

    def write_ccd(self, output_file=None, join=False):
        """Write the mmCIF file.

        Args:
            output_file (str, optional): The output file path. Defaults to None.
                If None, the output file will be {self.name}.cif
            join (bool, optional): Whether to join the output into a single line.
                Defaults to False.

        Returns:
            str: The path to the written file.
        """
        if self.mol is None:
            raise ValueError("Molecule not set. Cannot generate CIF file.")

        if output_file is None:
            if self.name is None:
                raise ValueError("Name not set and no output file specified.")
            output_file = f"{self.name}.cif"

        # Ensure all blocks are generated
        if self.comp is None or self.comp_atom is None or self.comp_bond is None:
            raise ValueError("CCD blocks not generated. Cannot write CIF file.")

        # Join the blocks into the final content
        content = f"{self.comp}\n{self.comp_atom}\n{self.comp_bond}"

        # Join into a single line if requested
        if join:
            content = content.replace("\n", "\\n")

        # Write the content to the output file
        with open(output_file, 'w') as f:
            f.write(content)

        return output_file

    def mol_to_ccd_comp_block(self, mol, comp_id):
        """Generate the composition block for the CIF file.

        Args:
            mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
            comp_id (str): Component ID for the CIF file

        Returns:
            str: The composition block.
        """
        # If molecule is not provided, return empty block
        if mol is None:
            return ""

        # Ensure comp_id is valid
        if comp_id is None:
            comp_id = "UNK"

        # Build comp block with single line entries instead of loop
        comp_lines = [
            "data_{0}".format(comp_id),
            "#",
            "_chem_comp.id {0}".format(comp_id),
            "_chem_comp.name \"{0}\"".format(self.iupac_name or comp_id),
            "_chem_comp.type non-polymer",
            "_chem_comp.formula \"{0}\"".format(self.formula or "Unknown"),
            "_chem_comp.mon_nstd_parent_comp_id ?",
            "_chem_comp.pdbx_synonyms ?",
            "_chem_comp.formula_weight {0:.3f}".format(self.weight or 0.0),
            "#"
        ]
        return "\n".join(comp_lines)

    def mol_to_ccd_comp_atom_block(self, mol, comp_id):
        """
        Convert RDKit molecule to custom CIF atom block format.

        Args:
            mol (rdkit.Chem.rdchem.Mol): RDKit molecule object with 3D coordinates
            comp_id (str): Component ID for the CIF file

        Returns:
            str: Content in custom CIF format
        """
        if mol is None:
            raise ValueError("Invalid RDKit molecule object")

        if mol.GetNumConformers() == 0:
            raise ValueError("Molecule must have 3D coordinates (at least one conformer)")

        # Initialize output list
        cif_lines = ["loop_",
                    "_chem_comp_atom.comp_id",
                    "_chem_comp_atom.atom_id",
                    "_chem_comp_atom.type_symbol",
                    "_chem_comp_atom.charge",
                    "_chem_comp_atom.pdbx_leaving_atom_flag",
                    "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
                    "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
                    "_chem_comp_atom.pdbx_model_Cartn_z_ideal"]

        # Track atom counters per element type
        element_counters = {}

        # Get conformer for 3D coordinates
        conf = mol.GetConformer()

        # Process atoms
        for atom in mol.GetAtoms():
            # Get atom properties
            idx = atom.GetIdx()
            element = atom.GetSymbol()
            charge = atom.GetFormalCharge()

            # Get 3D coordinates
            pos = conf.GetAtomPosition(idx)
            x, y, z = pos.x, pos.y, pos.z

            # Generate appropriate atom_id
            if element not in element_counters:
                element_counters[element] = 0
            element_counters[element] += 1

            if element == 'H':
                atom_id = f"H{element_counters[element]}"
            else:
                atom_id = f"{element}{element_counters[element]:02d}"

            # Add atom line to CIF content
            cif_lines.append(f"{comp_id} {atom_id} {element} {charge} N {x:.3f} {y:.3f} {z:.3f}")

        # Join all lines with newlines and return
        return "\n".join(cif_lines)

    def mol_to_ccd_comp_bond_block(self, mol, comp_id):
        """
        Convert RDKit molecule to custom CIF bond block format.

        Args:
            mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
            comp_id (str): Component ID for the CIF file

        Returns:
            str: Bond block in custom CIF format
        """
        if mol is None:
            raise ValueError("Invalid RDKit molecule object")

        # Generate atom IDs using the same convention as in mol_to_ccd_comp_atom_block
        atom_id_mapping = {}
        element_counters = {}

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            element = atom.GetSymbol()

            if element not in element_counters:
                element_counters[element] = 0
            element_counters[element] += 1

            if element == 'H':
                atom_id = f"H{element_counters[element]}"
            else:
                atom_id = f"{element}{element_counters[element]:02d}"

            atom_id_mapping[idx] = atom_id

        # Initialize bond block lines
        bond_lines = ["#",
                    "loop_",
                    "_chem_comp_bond.atom_id_1",
                    "_chem_comp_bond.atom_id_2",
                    "_chem_comp_bond.value_order",
                    "_chem_comp_bond.pdbx_aromatic_flag"]

        # Process bonds
        for bond in mol.GetBonds():
            # Get atom indices
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()

            # Get corresponding atom IDs
            begin_atom_id = atom_id_mapping[begin_atom_idx]
            end_atom_id = atom_id_mapping[end_atom_idx]

            # Determine bond type
            bond_type = bond.GetBondType()
            if bond_type == Chem.rdchem.BondType.SINGLE:
                bond_order = "SING"
            elif bond_type == Chem.rdchem.BondType.DOUBLE:
                bond_order = "DOUB"
            elif bond_type == Chem.rdchem.BondType.TRIPLE:
                bond_order = "TRIP"
            elif bond_type == Chem.rdchem.BondType.AROMATIC:
                bond_order = "SING" if not bond.GetIsAromatic() else "DOUB"
            else:
                bond_order = "SING"  # Default to single bond for other types

            # Determine aromaticity flag
            aromatic_flag = "Y" if bond.GetIsAromatic() else "N"

            # Add bond line to CIF content
            bond_lines.append(f"{begin_atom_id} {end_atom_id} {bond_order} {aromatic_flag}")

        # Add closing marker
        bond_lines.append("#")

        # Join all lines with newlines and return
        return "\n".join(bond_lines)


def generate_ccd_from_pdb(pdb_file, name=None) -> MyCCD:
    """Parse a PDB file and return a MyCCD object.

    Args:
        pdb_file (str): The path to the PDB file.
        name (str, optional): The name to use for the compound. Defaults to None.

    Returns:
        MyCCD: A MyCCD object.
    """
    # If name is not provided, use the file name without extension
    if name is None:
        name = os.path.splitext(os.path.basename(pdb_file))[0]

    # Read the PDB file content
    try:
        with open(pdb_file, 'r') as f:
            pdb_content = f.read()
        print(f"Converting PDB file: {pdb_file}")
    except Exception as e:
        raise ValueError(f"Could not read PDB file: {pdb_file}. Error: {str(e)}")

    # Convert PDB to RDKit molecule
    mol = Chem.MolFromPDBBlock(pdb_content)
    if mol is None:
        raise ValueError(f"Invalid PDB file: {pdb_file}")

    # Add Hs and minimize coords
    mol = Chem.AddHs(mol)
    Chem.MMFFOptimizeMolecule(mol)

    # Generate PDB block (we'll use RDKit here for consistency)
    pdb = Chem.MolToPDBBlock(mol)

    # Generate SMILES string
    smiles = Chem.MolToSmiles(mol)

    # Create a MyCCD object
    ccd = MyCCD(name, mol)

    # Ensure PDB and SMILES are set
    ccd.pdb = pdb
    ccd.smiles = smiles

    return ccd


def generate_ccd_from_smiles(smiles_input, name=None) -> MyCCD:
    """Parse a SMILES string and return a MyCCD object.

    Args:
        smiles_input (str): Either a SMILES string or path to a file containing a SMILES string.
        name (str, optional): The name to use for the compound. Defaults to None.
            If None and smiles_input is a file, the filename will be used.

    Returns:
        MyCCD: A MyCCD object.
    """
    smiles_string = ""

    # Check if input is a file or direct SMILES string
    if os.path.isfile(smiles_input):
        # Use filename as compound name if not specified
        name = name or os.path.splitext(os.path.basename(smiles_input))[0]

        try:
            with open(smiles_input, 'r') as f:
                smiles_string = f.read().strip()
            print(f"Converting SMILES file: {smiles_input}")
        except Exception as e:
            raise ValueError(f"Failed to read file: {smiles_input} - {e}")
    else:
        # Treat input as direct SMILES string
        smiles_string = smiles_input
        name = name or "MOL"
        print(f"Converting SMILES string: {smiles_string}")

    # Validate the SMILES
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles_string}")

    # Add Hs and generate coords
    mol = Chem.AddHs(mol)
    Chem.EmbedMolecule(mol)
    Chem.MMFFOptimizeMolecule(mol)

    # Generate PDB block
    pdb = Chem.MolToPDBBlock(mol)

    # Create a MyCCD object with the molecule
    ccd = MyCCD(name, mol)

    # Ensure PDB and SMILES are set
    ccd.pdb = pdb
    ccd.smiles = smiles_string

    return ccd