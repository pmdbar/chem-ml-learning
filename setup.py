from setuptools import setup, find_packages

setup(
    name="chemml",
    version="0.1.0",
    description="Cheminformatics utilities and CLI tools for small-molecule property computation.",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)

