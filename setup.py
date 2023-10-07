from distutils.core import setup
import setuptools

setup(
    name='msseq',
    version='1.0',
    author='Xiaohong Yuan',
    author_email='xyuan04@nyit.edu',
    packages = ['seq'],
    package_data = {
        'seq': ['statics/*.csv'],
    },
    test_suite="tests",
    entry_points={
        'console_scripts': [
            'seq = __main__:main'
        ]
    },
    install_requires = [
        'numpy>=1.16.2',
        'pandas>=1.3.5',
        'networkx>=2.2',
        'xlrd>=1.2.0',
        'matplotlib>=3.1.1',
        'statsmodels>=0.9.0',
        'scipy>=1.3.0',
        'loguru>=0.2.4'
    ],
    python_requires=">=3.5",
)
