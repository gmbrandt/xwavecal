from setuptools import setup, find_packages


with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='echelle',
      author=['G. Mirek Brandt, Curtis McCully'],
      version='0.1.0',
      python_requires='>=3.5',
      packages=find_packages(),
      package_dir={'echelle': 'echelle'},
      setup_requires=['pytest-runner'],
      install_requires=requirements,
      tests_require=['pytest>=3.5'],
      entry_points={'console_scripts': ['echelle_reduce_dir=echelle.main:run',
                                        'echelle_reduce=echelle.main:reduce_data']})
