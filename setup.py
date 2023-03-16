from setuptools import setup

from os import path
here = path.abspath(path.dirname(__file__))
def readme(file):
  with open(path.join(here, 'README.md')) as fh:
      long_description_text = fh.read()
  return(long_description_text)

if __name__ == "__main__":
	setup(name='pastrami', version='0.9.4',description='Pastrami is a novel, scalable computational algorithm for rapid human ancestry estimation at population-, subcontinental- and continental-levels. Pastrami works on two key methodologies: exact haplotype matching and non-negative least square (NNLS) optimization.',
python_requires='>=3.8',
zip_safe=False,
long_description=readme('README.md'),
long_description_content_type="text/markdown",
install_requires=['numpy','pathos','scipy','pandas', 'plink'], 
author = 'Jordan Lab',
  author_email = 'king.jordan@biology.gatech.edu',
  url = 'https://github.com/healthdisparities/pastrami',
keywords = ['ancestry', 'NNLS', 'haplotype matching'],
classifiers = [
      'Programming Language :: Python :: 3.8',
  ],
scripts = ['pastrami.py'],
	      py_modules = [],
   )
