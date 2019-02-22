import os
from setuptools import find_packages, setup


def extract_version():
    """
    Extracts version values from the main matplotlib __init__.py and
    returns them as a dictionary.
    """
    with open('NETSeq/__init__.py') as fd:
        for line in fd.readlines():
            if (line.startswith('__version__')):
                exec(line.strip())
    return locals()["__version__"]


def get_package_data():
    return {'NETSeq': ['examples/*.html', 'examples/*.txt', 'examples/*.ipynb']}

setup(
    name="NETSeq",
    version=extract_version(),
    author="Nudler Lab",
    author_email="ilya.shamovsky@gmail.com",
    url="https://github.com/NudlerLab/NET_analysis",
    license="MIT",
    packages=find_packages(),
    package_dir={"NETSeq": "NETSeq"},
    package_data=get_package_data(),
    description="NET-Seq data analysis pipeline",
    # run pandoc --from=markdown --to=rst --output=README.rst README.md
    long_description=open("README.rst").read(),
    # numpy is here to make installing easier... Needs to be at the
    # last position, as that's the first installed with
    # "python setup.py install"
    install_requires=["numpy",
                    "cutadapt",
                    "bash_kernel",
                      "pandas >= 0.16.0",
                      "seaborn",
                      "jupyter"],
    classifiers=['Intended Audience :: Science/Research',
                 'Programming Language :: Python',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Visualization',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: Unix',
                 'Operating System :: MacOS',
                 'Programming Language :: Python :: 2',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Programming Language :: Python :: 3.6'],
    zip_safe=False)
