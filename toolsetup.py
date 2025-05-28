import os
import subprocess
from typing import List
from enum import Enum
from abc import ABC, abstractmethod
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from pybind11.setup_helpers import Pybind11Extension, build_ext
import requests

EIGEN_DIRECTORY = os.path.abspath("src/core/math/eigen")

def install_eigen() -> None: 
    if not os.path.exists(EIGEN_DIRECTORY):
        os.makedirs(EIGEN_DIRECTORY, exist_ok=True)
        eigen_url = 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz'
        eigen_tar = 'eigen.tar.gz'
        # Use requests to download
        response = requests.get(eigen_url, stream=True)
        response.raise_for_status()
        with open(eigen_tar, 'wb') as f:
            f.write(response.content)
        subprocess.check_call(['tar', 'xzf', eigen_tar, '-C', EIGEN_DIRECTORY, '--strip-components=1'])
        os.remove(eigen_tar)

class CustomBuildExt(build_ext):
    def run(self):
        install_eigen()
        build_ext.run(self)
    
class ArbitrageBuild(ABC):

    def __init__(self) -> None:
        pass
    
    @abstractmethod
    def get_build_name(self) -> str: pass

    @abstractmethod
    def get_modules_names(self) -> List[str]: pass

    @abstractmethod
    def get_binding_file_name(self) -> str: pass 

    @abstractmethod
    def get_include_dirs(self) -> List[str]: pass

    @abstractmethod
    def get_library_dirs(self) -> List[str]: pass
    
    @abstractmethod
    def get_libraries(self) -> List[str]: pass 

    @abstractmethod
    def get_build_ext(self) -> build_ext: pass 

    @abstractmethod
    def get_setup_install_requires(self) -> List[str]: pass

    def get_package_name(self) -> str:
        return 'lightarbitrage'  
    
    def get_extra_compile_args(self) -> List[str]:
        return ["-std=c++17"]
    
    @staticmethod
    def list_files(module_name: str) -> List[str]:
        output = []
        endswith = '.cpp'
        directory = os.getcwd() + '/src/'+ module_name +'/'
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith(endswith):
                    output.append('src/'+os.path.relpath(os.path.join(root, file), 'src'))
        return output

    def __get_modules_names(self) -> List[str]:
        return self.get_modules_names() + ['core']

    def get_cpp_files(self) -> List[str]:
        output = []
        for module_name in self.__get_modules_names():
            output += self.list_files(module_name)
        return output
    
    def get_binding_file(self) -> str:
        #return os.getcwd() + '/src/pybindings/' + self.get_binding_file_name()
        return 'src/pybindings/' + self.get_binding_file_name()
    
    def __get_include_dirs(self) -> List[str]:
        return [EIGEN_DIRECTORY] + self.get_include_dirs()

    def get_extension(self) -> Pybind11Extension:
        files = self.get_cpp_files()
        files.append(self.get_binding_file())
        return Pybind11Extension(
                self.get_package_name(),
                sources=files,
                include_dirs = self.__get_include_dirs(), 
                extra_compile_args=self.get_extra_compile_args(),
                library_dirs=self.get_library_dirs(),
                libraries=self.get_libraries(),
            )
    
    def get_setup_kwargs(self) -> dict:
        return {
            'ext_modules' : [self.get_extension()],
            'cmdclass' : {'build_ext': self.get_build_ext()},
            'extras_require': {self.get_build_name() : self.get_setup_install_requires()},
        }

class CoreBuild(ArbitrageBuild):

    def get_build_name(self) -> str: return 'core'

    def get_modules_names(self) -> List[str]: return []

    def get_binding_file_name(self) -> str: return 'py_core_arbitrage.cpp'

    def get_include_dirs(self) -> List[str]: return []
    
    def get_library_dirs(self) -> List[str]: return []
    
    def get_libraries(self) -> List[str]: return []
    
    def get_build_ext(self) -> build_ext: return CustomBuildExt
    
    def get_setup_install_requires(self) -> List[str]: return []
    
build_list = [CoreBuild()]