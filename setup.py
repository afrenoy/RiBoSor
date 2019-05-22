import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="RiBoSor",
    url="https://github.com/afrenoy/RiBoSor",
    version="1.0",
    author="Antoine Frenoy",
    author_email="antoine.frenoy@pasteur.fr",
    description="Creates overlapping reading frames within an existing gene",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=["Programming Language :: Python :: 3", "License :: OSI Approved :: GNU General Public License v3 (GPLv3)", "Topic :: Scientific/Engineering :: Bio-Informatics"],
    install_requires=["biopython"],
    entry_points={'console_scripts': ["ribosor = ribosor.ribosor:main"]},
    packages=["ribosor"],
    include_package_data=True,
    )
