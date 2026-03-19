from setuptools import setup, find_packages

setup(
    name="Spiky",  # Change this if you want a different package name
    version="0.0.0",  # Start version, update as needed
    author="Your Name",
    author_email="your.email@example.com",
    description="A short description of the Spiky project",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/J35P312/Spiky",
    packages=find_packages(),  # Automatically find packages in the repo
    python_requires=">=3.8",  # Adjust based on your project
    install_requires=[
        "pysam>=0.20.0",
        "numpy>=1.21.0",
        "pandas>=1.4.0",
        "scipy>=1.8.0",
        "scikit-learn>=1.1.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Update if your license differs
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "spiky=spiky.main:main",  # Adjust if your main script/module differs
        ],
    },
)
