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
    version='1.0.0',
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
    python_requires='~=3.5.2'  # no commitment to Python 4 yet...
)

