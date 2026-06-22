from setuptools import setup, find_packages

setup(
    name="shisma",
    version="1.0.0",
    description="Shape-driven Inference of Significant celltype-specific Subnetworks",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "networkx",
        "statsmodels",
        "scikit-learn",
        "fasttreeshap",
        "scipy",
        "joblib"
    ],
    entry_points={
        "console_scripts": [
            "shisma=shisma.cli:main",
        ]
    },
)