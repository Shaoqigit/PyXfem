from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
  long_description = f.read()

setup(
    name='PyAcoustiX',
    version='0.9.13',
    description='A pure Python package/software for acoustic analysis',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Shaoqi WU',
    author_email='shaoqiwu@outlook.com',
    url='https://github.com/Shaoqigit/PyXfem',
    entry_points={
        'console_scripts': [
            'sacous = SAcouS.__main__:main',
        ],
    },
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'meshio',
    # Add any other dependencies here
    ],
    install_optional=[
        'pymls',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.8')
