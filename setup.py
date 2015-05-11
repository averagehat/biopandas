from setuptools import setup, find_packages

import your_project

setup(
    name = your_project.__projectname__,
    version = your_project.__release__,
    packages = find_packages(),
    author = your_project.__authors__,
    author_email = your_project.__authoremails__,
    description = your_project.__description__,
    license = "GPLv2",
    keywords = your_project.__keywords__,
    entry_points = {
        'console_scripts': [
            'your_project = your_project.your_project:main'
        ],
    },
)
