import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lsclib", # Replace with your own username
    version="0.1.2",
    author="Duncan E. Smith",
    author_email="smithd24@rpi.edu",
    description=("Library for evaluating the performance of LSCs"),
    long_description="""
    lsclib is a python-based repository hosted on GitHub that employs the 
    Monte Carlo ray-tracing method of radiative transport to effectively model
    LSCs. lsclib hopes to short-circuit the learning curve associated with
    breaking into the field, and present results in both an academic and 
    commercial context. This repository will continue to become more 
    sophisticated, but for now relies heavily upon the paper submitted for 
    publishing entitled "An open-source Monte Carlo ray-trace simulation tool 
    for luminescent solar concentrators with validation studies employing
    scattering phosphor films."
    """,
    long_description_content_type="text/markdown",
    url="https://github.com/duncanesmith/lsclib",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)