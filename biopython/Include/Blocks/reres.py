import os
import shutil

from HorusAPI import PluginBlock, VariableTypes, PluginVariable, VariableGroup

input_single_structure_id = "single_structure"
input_single_structure = VariableGroup(
    id=input_single_structure_id,
    name="Structure to reset",
    description="PDB file path to the structure to reset",
    variables=[
        PluginVariable(
            id=input_single_structure_id,
            name="PDB file",
            description="PDB file path to the structure to reset",
            type=VariableTypes.FILE,
            allowedValues=["pdb"],
        ),
    ],
)

folder_to_reset_id = "folder_to_reset"
folder_to_reset = VariableGroup(
    id=folder_to_reset_id,
    name="Folder to reset",
    description="Folder containing PDBs to be reset to the reference",
    variables=[
        PluginVariable(
            id=folder_to_reset_id,
            name="Folder to reset",
            description="Folder containing PDBs to be reset to the reference",
            type=VariableTypes.FOLDER,
        ),
    ],
)


starting_res_var = PluginVariable(
    id="starting_res_var",
    name="Starting residue",
    description="The index of the starting residue",
    type=VariableTypes.NUMBER,
    defaultValue=1,
)

output_pdb = PluginVariable(
    id="pdb_path",
    name="Reset PDBs",
    description="The folder containing the reset PDB files",
    type=VariableTypes.FOLDER,
)



def reset_residue_numbers(pdb_file, start_number=1, folder: str = "out"):
    from Bio import PDB
    # Initialize the PDB parser
    parser = PDB.PPBuilder()
    structure = PDB.PDBParser(QUIET=True).get_structure('structure', pdb_file)

    # Start numbering residues from start_number
    current_residue_number = start_number

    # Iterate through each model in the structure (usually one model in PDB)
    for model in structure:
        for chain in model:
            for residue in chain:
                # Reset the residue number
                residue.id = (' ', current_residue_number, ' ')
                current_residue_number += 1

    # Write the modified structure to a new PDB file
    io = PDB.PDBIO()
    io.set_structure(structure)

    output_filename = os.path.basename(pdb_file).split(".")[0] + "_reset.pdb"

    io.save(os.path.join(folder, output_filename))

    return output_filename


def reset(block: PluginBlock):

    start_number = block.variables[starting_res_var.id]

    pdbs_to_reset = []
    if block.selectedInputGroup == input_single_structure_id:
        pdbs_to_reset = [block.inputs[input_single_structure_id]]
    else:
        folder = block.inputs[folder_to_reset_id]
        for path in os.listdir(folder):
            pdbs_to_reset.append(os.path.join(folder, path))

    # Make the reset folder
    output_folder = "reset_pdb_folder"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    os.makedirs(output_folder, exist_ok=True)

    for pdb in pdbs_to_reset:
        output_file_name_value = reset_residue_numbers(pdb, start_number=start_number, folder=output_folder)


        print(f" - Reset PDB saved as: {output_file_name_value}")

    block.setOutput(output_pdb.id, output_folder)


reset_block = PluginBlock(
    id="reset",
    name="Reset Residues in PDBs",
    description="Reset residue numbers to start at the given number",
    action=reset,
    inputGroups=[input_single_structure, folder_to_reset],
    variables=[starting_res_var],
    outputs=[output_pdb],
)
