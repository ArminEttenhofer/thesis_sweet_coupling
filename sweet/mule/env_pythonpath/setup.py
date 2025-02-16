#! /usr/bin/env python3

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

print(setuptools.find_packages())

setuptools.setup(
    name="mule",
    version="0.0.1",
    author="Martin Schreiber et al.",
    author_email="schreiberx@gmail.com",
    description="MULE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['mule'],
    python_requires='>=3.6',
)

