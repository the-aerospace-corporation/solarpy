import setuptools
from setuptools import setup

setup(
        name='solarpy',
        version='0.1a',
        packages=setuptools.find_packages(),
        package_data={'solarpy': ['data/*.txt', 'data/*.xlsx']},
        url='https://aerosource2.aero.org/bitbucket/scm/~dw28596/solarpy.git',
        license='GNU GPLv3',
        author='Don Walker',
        author_email='don.walker@aero.org',
        description='Data analysis tools for solar cells',
        install_requires=['numpy',
                          'scipy',
                          'pandas',
                          'matplotlib',
                          'julian',
                          'profilehooks',
                          'openpyxl'],

)
