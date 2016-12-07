import os
from numpy.distutils.core import setup, Extension

utils    = Extension(name = 'OMG.utils',
                             sources = ['OMG/f90/utils.f90'])


setup(
    name = "OMG",
    version = "1.0",
    author = "Raphael Dussin",
    author_email = "raphael.dussin@gmail.com",
    description = ("Ocean Model Google earth plots in python " ),
    license = "GPLv3",
    keywords = "ocean modeling, visualization",
    url = "",
    packages=['OMG'], 
    ext_modules = [utils]
)


