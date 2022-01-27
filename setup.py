#!/usr/bin/env python
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GELATO",
    version="2.5",
    author="Raphael Hviding",
    author_email="raphael.hviding@gmail.com",
    description="Galaxy/AGN Emission Line Analysis TOol ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TheSkyentist/GELATO",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3',
)
