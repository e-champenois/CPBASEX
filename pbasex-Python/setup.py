from setuptools import setup

if __name__ == '__main__':
	setup(
		  name='pbasex',
		  keywords='pbasex, Abel, inversion, cpbasex',
		  packages=['pbasex'],
		  setup_requires=['cython'],
		  install_requires=['numpy', 'h5py', 'quadrant'],
		  version='1.2.2',
		  description='pBASEX algorithm without polar rebinning.',
		  url='https://github.com/e-champenois/CPBASEX',
		  author='Elio Champenois',
		  author_email='elio.champenois@gmail.com',
		  license='MIT',
		  zip_safe=False
		  )