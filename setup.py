from setuptools import setup, find_packages

setup(
    name='Proj_VCB',  # Name of your package
    version='0.1.0',  # Initial release version
    description='A package for analyzing, filtering and visualizing cosmetic creams with RDKit. and panda',
    long_description=open('README.md').read(),  # Read the long description from your README file
    long_description_content_type='text/markdown',  # Specify the content type of the long description
    url='https://github.com/victoirelang/Project',  # URL to your project repository
    author='Victoire Lang, Blanche Billarant, Clara Moret',
    author_email='victoire.lang@epfl.ch, blanche.billarant-laubignac@epfl.ch, clara.moret@epfl.ch',
    license='MIT',  # Choose your license, e.g., MIT
    packages=find_packages(where='src'),  # Include all packages under the src directory
    package_dir={'': 'src'},  # Tell setuptools to look in the src directory for packages
    include_package_data=True,  # Include package data specified in MANIFEST.in
    install_requires=[
        'pandas',  # List all dependencies
        'rdkit',
        'IPython',
    ],
    extras_require={
        'dev': [
            'pytest',
            'sphinx',
            'flake8',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',  # Define the development status
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.8',  # Specify the Python versions you support
)
