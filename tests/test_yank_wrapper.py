import unittest
from unittest.mock import Mock

from cwl_adapters.yank import yank_wrapper


class TestYankWrapper(unittest.TestCase):

    def test_check_options(self):
        options = {}
        args = Mock()
        yank_wrapper.check_options(options, args)
        self.assertEqual(options, {
            'verbose': True,
            'resume_setup': True,
            'resume_simulation': True,
            'start_from_trailblaze_samples': False,
            'checkpoint_interval': 10,
            'switch_phase_interval': 10,
            'number_of_equilibration_iterations': 25,
            'default_number_of_iterations': 10
        })

    def test_make_molecules_tag(self):
        args = Mock()
        args.input_receptor_path = 'receptor_path'
        args.input_ligand_path = 'ligand_path'
        result = yank_wrapper.make_molecules_tag(args)
        self.assertEqual(result, {
            'receptor_name': {'filepath': 'receptor_path'},
            'ligand_name': {'filepath': 'ligand_path', 'antechamber': {'charge_method': 'bcc'}}
        })

    def test_make_systems_tag_group1(self):
        result = yank_wrapper.make_systems_tag_group1()
        self.assertEqual(result, {
            'systems': {
                'system_name': {
                    'receptor': 'receptor_name',
                    'ligand': 'ligand_name',
                    'solvent': 'spce_50mM_acetate_55',
                    'leap': {'parameters': ['oldff/leaprc.ff99SBildn', 'leaprc.gaff']}
                }
            }
        })

    def test_make_systems_tag_group2(self):
        args = Mock()
        args.input_complex_crd_path = 'complex_crd_path'
        args.input_ligand_crd_path = 'ligand_crd_path'
        result = yank_wrapper.make_systems_tag_group2('complex_top_filename', 'ligand_top_filename', args)
        self.assertEqual(result, {
            'systems': {
                'system_name': {
                    'phase1_path': ['complex_top_filename', 'complex_crd_path'],
                    'phase2_path': ['ligand_top_filename', 'ligand_crd_path'],
                    'ligand_dsl': 'resname MOL',
                    'solvent': 'spce_50mM_acetate_55',
                    'gromacs_include_dir': '/miniconda/share/gromacs/top/'
                }
            }
        })

    def test_make_experiment_tag(self):
        result = yank_wrapper.make_experiment_tag()
        self.assertEqual(result, {
            'system': 'system_name',
            'sampler': 'repex',
            'options': {
                'temperature': '302.15*kelvin'
            },
            'protocol': 'binding-auto'
        })

    def test_unzip_topology_files(self):
        with self.assertRaises(FileNotFoundError):
            yank_wrapper.unzip_topology_files('nonexistent.zip')

    def test_cli(self):
        parser = yank_wrapper.cli()
        self.assertEqual(type(parser), argparse.ArgumentParser)

if __name__ == '__main__':
    unittest.main()
