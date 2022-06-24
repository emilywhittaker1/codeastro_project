from os import mkdir
from setuptools import setup, find_packages
import os

plotdir = os.path.join(os.path.dirname(__file__), "Plot")

if os.path.exists(plotdir):
  os.mkdir(plotdir)

setup(
name="differential_photometry",
version="0.2",
packages=find_packages()
)


