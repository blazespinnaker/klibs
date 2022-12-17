import io
from setuptools import setup, find_packages
setup(
    name='klib',
    version='1.7.0',
    description='utils for k',
    author='blaze',
    author_email='blaze@gmail.com',
    url='https://github.com/blazespinnaker/klibs',
    license='MIT',
    install_requires=['numpy','pandas', 'pandas_profiling'],
    py_modules=['klibs'],
    packages=find_packages(),
    long_description_content_type='text/markdown',
    long_description=io.open('README.md', encoding='utf-8').read()
)
