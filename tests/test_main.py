import unittest

from examples.scripts.python_cwl_driver import import_python_file


class TestMain(unittest.TestCase):

    def test_import_python_file(self):
        # Test if the function correctly imports a python file
        module = import_python_file('workflow_types', 'workflow_types.py')
        self.assertIsInstance(module, ModuleType)

    def test_import_python_file_invalid(self):
        # Test if the function raises an error when trying to import a non-existent file
        with self.assertRaises(ModuleNotFoundError):
            import_python_file('non_existent_module', 'non_existent_file.py')

if __name__ == '__main__':
    unittest.main()
