"""
    setup.py
    Dylan H. Ross
    2020/01/17

    description:
        setup script for lipydomics package
"""


import setuptools


# use the README file contents for long description
with open('README.md', 'r') as f:
    long_desc = f.read()


setuptools.setup(
    name='lipydomics',
    version='1.6.8',
    author='Dylan H. Ross',
    author_email='dhross92@uw.edu',
    description='a library for streamlining lipidomics data analysis',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://github.com/dylanhross/lipydomics',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    include_package_data=True,  # include any data files within the package when building the distribution
    python_requires='>=3.8',  # support Python3.8 and up
    install_requires=[  # install or upgrade dependencies
        'matplotlib>=3.1.3',
        'numpy>=1.18.1',
        'pandas>=1.0.1',
        'scikit-learn>=0.24.1',
        'scipy>=1.4.1'
    ]
)

