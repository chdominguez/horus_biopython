import os
import shutil

from HorusAPI import PluginBlock, VariableTypes, PluginVariable, VariableGroup

input_structure_ref = PluginVariable(
    id="structure_1",
    name="Reference structure",
    description="The referebce structure to align the others to",
    type=VariableTypes.FILE,
    allowedValues=["pdb"],
)

input_single_structure_id = "single_structure"
input_single_structure = VariableGroup(
    id=input_single_structure_id,
    name="Structure to align",
    description="PDB file path to the structure to align",
    variables=[
        input_structure_ref,
        PluginVariable(
            id=input_single_structure_id,
            name="Second structure",
            description="The second structure to align",
            type=VariableTypes.FILE,
            allowedValues=["pdb"],
        ),
    ],
)

folder_to_align_id = "folder_to_align"
folder_to_align = VariableGroup(
    id=folder_to_align_id,
    name="Folder to align",
    description="Folder containing PDBs to be aligned to the reference",
    variables=[
        input_structure_ref,
        PluginVariable(
            id=folder_to_align_id,
            name="Folder to align",
            description="Folder containing PDBs to be aligned to the reference",
            type=VariableTypes.FOLDER,
        ),
    ],
)

output_file_name = PluginVariable(
    id="filename",
    name="Output file suffix",
    description="The suffix of the aligned PDB files",
    type=VariableTypes.STRING,
    defaultValue="_aligned",
)

output_pdb = PluginVariable(
    id="pdb_path",
    name="Aligned PDBs",
    description="The folder containing the aligned PDB files",
    type=VariableTypes.FOLDER,
)


def use_biopython_to_align(ref_pdb, mov_pdb, filename):

    # https://gist.github.com/andersx/6354971

    # The MIT License
    #
    # Copyright (c) 2010-2016 Anders S. Christensen
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    # THE SOFTWARE.

    import Bio.PDB

    # Select what residues numbers you wish to align
    # and put them in a list
    start_id = 1
    end_id = 70
    atoms_to_be_aligned = range(start_id, end_id + 1)

    # Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    # Get the structures
    ref_structure = pdb_parser.get_structure("reference", ref_pdb)
    sample_structure = pdb_parser.get_structure("samle", mov_pdb)

    # Use the first model in the pdb-files for alignment
    # Change the number 0 if you want to align to another structure
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]

    # Make a list of the atoms (in the structures) you wish to align.
    # In this case we use CA atoms whose index is in the specified range
    ref_atoms = []
    sample_atoms = []

    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
        for ref_res in ref_chain:
            # Check if residue number ( .get_id() ) is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                # Append CA atom to list
                ref_atoms.append(ref_res["CA"])

    # Do the same for the sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_atoms.append(sample_res["CA"])

    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # Print RMSD:
    print(f" - RMSD: {super_imposer.rms}")

    # Save the aligned version of 1UBQ.pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure)
    io.save(filename)


def align(block: PluginBlock):
    ref_pdb = block.inputs[input_structure_ref.id]

    pdbs_to_align = []
    if block.selectedInputGroup == input_single_structure_id:
        pdbs_to_align = [block.inputs[input_single_structure_id]]
    else:
        folder = block.inputs[folder_to_align_id]
        for path in os.listdir(folder):
            pdbs_to_align.append(os.path.join(folder, path))

    # Make the alignment folder
    output_folder = "aligned_folder"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    os.makedirs(output_folder, exist_ok=True)

    for pdb in pdbs_to_align:
        suffix = block.variables[output_file_name.id] or "_aligned"
        output_file_name_value = os.path.basename(pdb).split(".")[0] + f"{suffix}.pdb"
        output_file_name_value = os.path.join(output_folder, output_file_name_value)

        use_biopython_to_align(ref_pdb, pdb, output_file_name_value)

        print(f" - Aligned PDB saved as: {output_file_name_value}")

    block.setOutput(output_pdb.id, output_folder)


align_block = PluginBlock(
    id="align",
    name="Align PDBs",
    description="Align 3D PDB structures and print their RMSD alignment",
    action=align,
    inputGroups=[input_single_structure, folder_to_align],
    variables=[output_file_name],
    outputs=[output_pdb],
)
