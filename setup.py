import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sp-pharkin",  # Replace with your own username
    version="0.0.3",
    author="Stavros Pitoglou",
    author_email="spitoglou@biomed.ntua.gr",
    description="Pharmacokinetics Helper Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/spitoglou/sp-pharkin",
    packages=setuptools.find_packages(),
    install_requires=[
        'pint',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
