from setuptools import setup, find_packages

setup(
    name='Proj_VCB',  
    version='0.1.0',  
    description='A package for analyzing, filtering and visualizing cosmetic creams with RDKit. and panda',
    long_description=open('README.md').read(),  
    long_description_content_type='text/markdown',  
    url='https://github.com/victoirelang/Project',  
    author='Victoire Lang, Blanche Billarant, Clara Moret',
    author_email='victoire.lang@epfl.ch, blanche.billarant-laubignat@epfl.ch, clara.moret@epfl.ch',
    license='MIT',  
    packages=find_packages(where='src'),  
    package_dir={'': 'src'},  
    include_package_data=True,  
    install_requires=[
        'pandas',  
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
        'Development Status :: 3 - Alpha',  
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.8',  
)
