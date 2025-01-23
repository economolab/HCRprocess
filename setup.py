import setuptools
from setuptools import setup

docs_deps = [
    'sphinx>=3.0',
    'sphinxcontrib-apidoc',
    'sphinx_rtd_theme',
    'sphinx-argparse',
]

setup(
    name="HCRprocess", license="BSD", author="Jack Vincent",
    author_email="jackv@bu.edu",
    description="HCR-mFISH processing pipeline", long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/economolab/HCRprocess", setup_requires=[
        'pytest-runner',
        'setuptools_scm',
    ], packages=setuptools.find_packages(), use_scm_version=True,
    install_requires=install_deps, tests_require=['pytest'], extras_require={
        'docs': docs_deps,
    }, include_package_data=True, classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ), entry_points={'console_scripts': ['cellpose = cellpose.__main__:main']})
