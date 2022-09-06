from distutils.core import setup


# Build the f2py fortran extension
# --------------------------------
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

#import os
#os.chdir('Fortran/bhmie/')
#os.system('gfortran -c -O3 -Wall -Werror -fPIC bhmie.f')
#os.system('gcc -shared -o libbhmie.so bhmie.o')
#os.chdir('../../')
#print("here is my setup thing")

flib = Extension(name = 'retLib.bhmie',
                 extra_compile_args = ['-O3'],
                 sources = ['Fortran/bhmie/mieRoutines.f90'],
                 include_dirs=[],
                 libraries=['bhmie'],
                 library_dirs=['/home/grecu/packaging_dir/Fortran/bhmie/'],
                 )
calg = Extension(name = 'retLib.cAlgLite',
                 extra_compile_args = ['-O3','-fopenmp'],
                 extra_f90_compile_args = ['-O3','-fopenmp'],
                 extra_f77_compile_args = ['-O3','-fopenmp'],
                 sources = ['/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/readTables_new.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/readTables_nonsph.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/bisection.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/hbprof_new.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/absorption3D.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/rosen.f',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/emissivity-sp.f',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/gcloud.f',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/tbCalc.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/radtran_py_wrap.f90',\
                            '/home/grecu/packaging_dir/Fortran/cAlgLite/src_f90/dfr_retr.f90'
                            
                 ],
                 include_dirs=[],
                 libraries=['rte','gomp'],
                 library_dirs=['/home/grecu/packaging_dir/Fortran/cAlgLite/'],
                 )
setup(
    name = 'retLib',
    description       = "example package using fortran",
    author            = "Mircea Grecu",
    packages = ["retLib"],
    ext_modules = [flib,calg],
    package_data = {
        'retLib': [],
        'retLib.bhmie': [],
        'retLib': ['supportData/nw_dm_dsd/*', \
                   'supportData/TablesN/*', 'supportData/AncData/*',\
                   'supportData/zku_dm.txt',\
                   'librte.so',
                   'libbhmie.so']
        },
    )
