from setuptools import setup, find_packages

import python_template

setup(
    name = python_template.__projectname__,
    version = python_template.__release__,
    packages = find_packages(),
    author = python_template.__authors__,
    author_email = python_template.__authoremails__,
    description = python_template.__description__,
    license = "GPLv2",
    keywords = python_template.__keywords__,
    entry_points = {
        'console_scripts': [
            'python_template = python_template.python_template:main'
        ],
    },
)
