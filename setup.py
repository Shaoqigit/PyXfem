from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name='PyAcoustiX',
    version='0.9.9',
    description='A Python package for acoustic analysis',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Shaoqi WU',
    author_email='shaoqiwu@outlook.com',
    url='https://github.com/Shaoqigit/PyXfem',
    entry_points={
        'console_scripts': [
            'sacous = SAcouS.main:main',
        ],
    },
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'numba',
        # Add any other dependencies here
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    python_requires='>=3.8'
)