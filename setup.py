from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

pypi_classifiers = [
    "Programming Language :: Python",
    "Development Status :: 3 - Alpha",
    "Environment :: Console",
    "Operating System :: OS Independent",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    "pandas>=0.20.3",
    'biopython>=1.70',
]

desc = """
Map TIR-pHMM models to genomic sequences for annotation of MITES and complete DNA-Transposons.
"""

setup(name='tirmite',
      version='0.1.0',
      description=desc,
      long_description=readme(),
      url='https://github.com/Adamtaranto/TIRmite',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['tirmite'],
      classifiers=pypi_classifiers,
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'tirmite=tirmite.command_line:main',
        ],
    },
    )
