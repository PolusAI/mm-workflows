import unittest

import docker


class DockerfileSminaTest(unittest.TestCase):
    def setUp(self):
        self.client = docker.from_env()

    def test_build_image(self):
        image, _ = self.client.images.build(path="examples/scripts/", dockerfile="Dockerfile_smina")
        self.assertIsNotNone(image)

    def test_packages_installed(self):
        image, _ = self.client.images.build(path="examples/scripts/", dockerfile="Dockerfile_smina")
        container = self.client.containers.run(image, command="bash -c 'dpkg -l | grep -E \"git|build-essential|libboost-all-dev|libopenbabel-dev|libeigen3-dev|cmake\"'", remove=True)
        self.assertIn("git", container.decode())
        self.assertIn("build-essential", container.decode())
        self.assertIn("libboost-all-dev", container.decode())
        self.assertIn("libopenbabel-dev", container.decode())
        self.assertIn("libeigen3-dev", container.decode())
        self.assertIn("cmake", container.decode())

    def test_smina_compiled(self):
        image, _ = self.client.images.build(path="examples/scripts/", dockerfile="Dockerfile_smina")
        container = self.client.containers.run(image, command="bash -c '/smina-code/build/smina --help'", remove=True)
        self.assertIn("--help                        display usage summary", container.decode())

if __name__ == "__main__":
    unittest.main()
