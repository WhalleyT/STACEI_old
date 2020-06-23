import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="STACEI-WhalleyT",
    version="0.1",
    author="Tom Whalley",
    author_email="twhalley93@gmail.com",
    description="Python package for STACEI: structural analysis for TCR-pMHC complexes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/whalleyt/stacei",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU License",
        "Operating System :: Linux",
    ],
    entry_points={"console_scripts": [
        "STACEI = stacei.STACEI:main"
    ]},
    python_requires='>=3.6',
    install_requires=open("requirements.txt").read(),
    dependency_links=[
        "https://github.com/oxpig/ANARCI#egg=anarci",
        "https://github.com/schrodinger/pymol-open-source#egg=pymol"
    ],
    package_data={'stacei': ['bin/web/*', 'bin/data/*', 'bin/executables/*', 'bin/R/*']}
)