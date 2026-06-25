from setuptools import setup, find_packages

setup(
    name="shisma",
    version="1.0.0",
    description="Subnetwork Analysis via Bag-of-Receptive-Fields (BORF) and SHAP",
    author="Bioinformatics Engineer",
    packages=find_packages(),
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "shisma=shisma.cli:main",
        ],
    },
)