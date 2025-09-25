from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="thermompnn",
    version="1.0.0",
    author="Henry Dieckhaus, Michael Brocidiacono, Nicholas Z. Randolph, Brian Kuhlman",
    author_email="kuhlman@email.unc.edu",
    description="A graph neural network for predicting changes in protein stability",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Kuhlman-Lab/ThermoMPNN",
    packages=find_packages(include=["thermompnn", "thermompnn.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=[
        "torch>=2.0.0",
        "torchvision",
        "torchaudio",
        "pytorch-lightning",
        "biopython",
        "wandb",
        "tqdm",
        "pandas",
        "numpy",
        "omegaconf",
        "joblib",
    ],
    include_package_data=True,
    package_data={
        "thermompnn": [
            "../models/*.pt",
            "../config.yaml",
            "../local.yaml",
            "../dataset_splits/*.pkl",
        ],
    },
)
