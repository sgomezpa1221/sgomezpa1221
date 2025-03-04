from setuptools import setup, find_packages

setup(
    name="InfluenzaAnalysis",
    version="1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython",
        "pandas",
        "seaborn",
        "matplotlib",
        "numpy"
    ],
    entry_points={
        "console_scripts": [
            "influenza-analysis=src.main:main"
        ]
    },
    author="Jonatha Agudelo & Sergio Gómez",
    description="Un paquete para el análisis del virus de la influenza.",
    license="MIT",
)
