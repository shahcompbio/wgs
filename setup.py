from setuptools import setup, find_packages
import versioneer


setup(
    name='wgs',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='variant_calling',
    author='Diljot Grewal',
    author_email='diljot.grewal@gmail.com',
    entry_points={'console_scripts': ['wgs = wgs.run:main', 'wgs_bamtofastq = wgs.workflows.realignment.bamtofastq:main']},
    package_data={'':['scripts/*.py', 'scripts/*.R', 'scripts/*.npz', 'scripts/*.pl', 'scripts/*.Rmd', "config/*.yaml", "data/*"]}
)
