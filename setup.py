from setuptools import setup, find_packages

setup(name='xwavecal',
      author='G. Mirek Brandt, Curtis McCully, Timothy Brandt',
      author_email='gmbrandt@ucsb.edu',
      version='0.1.11',
      python_requires='>=3.6',
      url='https://github.com/gmbrandt/xwavecal',
      packages=find_packages(),
      package_dir={'xwavecal': 'xwavecal'},
      package_data={'xwavecal': ['data/.*']},
      setup_requires=['pytest-runner'],
      install_requires=[
                        'astropy<=4.0',
                        'scipy<=1.3.2',
                        'sep<=1.10.0',
                        'numpy>=1.16',
      ],
      entry_points={'console_scripts': ['xwavecal_reduce_dir=xwavecal.main:run',
                                        'xwavecal_reduce=xwavecal.main:reduce_data']})
