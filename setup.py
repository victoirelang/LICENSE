from setuptools import setup, find_packages

setup(
    name="my_package",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # Ajoutez les dépendances ici, par exemple :
        # "numpy",
        # "pandas",
    ],
    author="Victoire, Clara, Blanche",
    author_email="victoire.lang@epfl.ch, blanche.billarant-laubignac@epfl.ch,clara.moret@epfl.ch",
    description="Une brève description de votre package",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/victoirelang/Project",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
