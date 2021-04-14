import setuptools

setuptools.setup(
    name="fragmenter",
    version="0.0.1",
    author="Robert R. Puccinelli",
    author_email="robert.puccinelli@outlook.com",
    description="Utilities for plasmid fragment generation.",
    url="https://github.com/Lab-Unofficial/Fragmenter",
    packages=setuptools.find_packages(exclude=["*.tests", "*.tests.*",
                                               "tests.*", "tests"]),
    install_requires=[
    ],
    test_suite="tests",
    classifiers=[
        "UCSF :: DeRisi",
    ],
)