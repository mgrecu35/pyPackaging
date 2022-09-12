from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor

from netCDF4 import Dataset
import numpy as np

def train_knn(nn):
    fname='Chase_et_al_2021_NN-master/Unrimed_simulation_wholespecturm_train_V2.nc'
    fh=Dataset(fname)
    zKu_train=fh["Z"][:]
    zKa_train=fh["Z2"][:]
    iwc_train=fh["IWC"][:]
    temp_train=fh["T_env"][:]
    Nw_train=fh["Nw"][:]
    Dm_train=fh["Dm"][:]
    n=zKu_train.shape[0]
    zKu_train+=np.random.randn(n)*3
    zKa_train+=np.random.randn(n)*3
    zL_train=np.array([zKu_train,zKa_train]).T
    
    knn=KNeighborsRegressor(n_neighbors=nn,weights='distance')
    knn_dm=KNeighborsRegressor(n_neighbors=nn,weights='distance')
    knn.fit(zL_train.data,iwc_train.data)
    knn_dm.fit(zL_train.data,1e3*Dm_train.data)
    return knn,knn_dm,zL_train,iwc_train,Nw_train

def train_dtree(nn):
    fname='Chase_et_al_2021_NN-master/Unrimed_simulation_wholespecturm_train_V2.nc'
    fh=Dataset(fname)
    zKu_train=fh["Z"][:]
    zKa_train=fh["Z2"][:]
    iwc_train=fh["IWC"][:]
    temp_train=fh["T_env"][:]
    Nw_train=fh["Nw"][:]
    Dm_train=fh["Dm"][:]
    n=zKu_train.shape[0]
    zKu_train+=np.random.randn(n)*3
    zKa_train+=np.random.randn(n)*3
    zL_train=np.array([zKu_train,zKa_train]).T
    
    dtree=DecisionTreeRegressor(max_depth=5)
    dtree_dm=DecisionTreeRegressor(max_depth=5)
    dtree.fit(zL_train.data,iwc_train.data)
    dtree_dm.fit(zL_train.data,1e3*Dm_train.data)
    return dtree,dtree_dm,zL_train,iwc_train,Nw_train





def get_validation():
    fname='Chase_et_al_2021_NN-master/Unrimed_simulation_wholespecturm_test_V2.nc'
    fh=Dataset(fname)
    zKu_train=fh["Z"][:]
    zKa_train=fh["Z2"][:]
    iwc_train=fh["IWC"][:]
    temp_train=fh["T_env"][:]
    Nw_train=fh["Nw"][:]
    Dm_train=fh["Dm"][:]
    n=zKu_train.shape[0]
    zKu_train+=np.random.randn(n)*3
    zKa_train+=np.random.randn(n)*3
    zL_train=np.array([zKu_train,zKa_train]).T
    return zL_train,iwc_train,Dm_train, temp_train
