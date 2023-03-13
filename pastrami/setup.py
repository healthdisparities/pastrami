from setuptools import setup

if __name__ == "__main__":
	setup(name='pastrami', version='0.9.0',description='Pastrami is a novel, scalable computational algorithm for rapid human ancestry estimation at population-, subcontinental- and continental-levels. Pastrami works on two key methodologies: exact haplotype matching and non-negative least square (NNLS) optimization.',python_requires='>=3.8',zip_safe=False,install_requires=['numpy','pathos','scipy','pandas', 'plink'])
