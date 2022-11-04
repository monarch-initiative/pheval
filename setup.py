from setuptools import setup, find_packages

setup(
    name="pheval",
    version="0.1.0",
    packages=["pheval"],
    entry_points={"console_scripts": ["pheval = pheval.cli:main"]},
)
