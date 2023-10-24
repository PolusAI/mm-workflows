import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_pdb', type=str)
parser.add_argument('--output_pdb', type=str)
args = parser.parse_args()


def clean_smina_pdb(input_pdb: str, output_pdb: str) -> None:
    """Clean smina pdb file (clean the pdb and rename the resname to LIG)

    Args:
        input_pdb (str): input pdb file
        output_pdb (str): output pdb file
    """
    with open(input_pdb, mode='r', encoding='utf-8') as f1, open(output_pdb, mode='w', encoding='utf-8') as f2:
        for line in f1.readlines():
            if line.startswith(('ATOM', 'HETATM', 'CONNECT', 'TER')):
                # https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf
                # residue name is present in column 18-20, since the column starts from 0, should be 17-19.
                # onionnet-v1 based on residue name (LIG) to recognize ligand structure.
                if len(line) >= 21:
                    f2.write(line[:17]+'LIG'+line[20:])
                else:
                    f2.write(line)


if __name__ == '__main__':
    clean_smina_pdb(args.input_pdb, args.output_pdb)
