from abc import ABCMeta, abstractmethod
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.decomposition import FactorAnalysis as sk_fa
import h5py
from ZIFA import ZIFA, block_ZIFA
import numpy as np
import pandas as pd

class UnsupervisedMethods(metaclass=ABCMeta):

    def __init__(self):
        self.matrix = None
        self.barcodes = None
        self.results = None

    def load_from_csv(self, input_file):
        pass

    def load_from_hdf5(self, input_file):
        h5f = h5py.File(input_file, 'r')
        self.matrix = h5f['matrix'].value
        self.barcodes = h5f['cell_attrs']['cell_names'].value
        h5f.close()
 

    def write_csv(self, output_file):
        df = pd.DataFrame(self.barcodes)
        df = pd.concat([df, pd.DataFrame(self.results)], axis=1)
        df.to_csv(output_file,header=False, index=False)

    @abstractmethod
    def apply(self):
        pass

    def log_normalize(self):
        self.matrix = np.log1p(self.matrix)

class Pca(UnsupervisedMethods):

    def __init__(self, n_components):
        UnsupervisedMethods.__init__(self)
        self.n_components = n_components

    def apply(self):
        pca = PCA(n_components = self.n_components)
        self.results = pca.fit_transform(self.matrix)

class Tsne(UnsupervisedMethods):

    def __init__(self, n_components, init):
        UnsupervisedMethods.__init__(self)
        self.n_components = n_components
        self.init = init

    def apply(self):
        tsne = TSNE(init=self.init, n_components=self.n_components)
        self.results = tsne.fit_transform(self.matrix)

class FactorAnalysis(UnsupervisedMethods):

    def __init__(self, n_components):
        UnsupervisedMethods.__init__(self)
        self.n_components = n_components

    def apply(self):
        fa = sk_fa(n_components = self.n_components)
        self.results = fa.fit_transform(self.matrix)

class Zifa(UnsupervisedMethods):

    def __init__(self, n_components, n_blocks):
        UnsupervisedMethods.__init__(self)
        self.n_components = n_components
        self.n_blocks = n_blocks

    def apply(self):
        Zhat, params = block_ZIFA.fitModel(self.matrix, n_components= self.n_components, n_blocks = self.n_blocks)
        self.results = Zhat

class Simlr(UnsupervisedMethods):
    
    def __init__(self,n_components, pca_components, n_neighbours, max_iter, threads):
        UnsupervisedMethods.__init__(self)
        self.n_components = n_components
        self.pca_components = pca_components
        self.n_neighbours = n_neighbours
        self.max_iter = max_iter
        self.threads = threads

    def apply(self):
        X = self.matrix
        # Selecting 500 features with PCA
        if X.shape[1]>500:
        # fast_pca assumes the number of cells > 500 therefore try-catch
            try:
                X = SIMLR.helper.fast_pca(X,self.pca_components)
            except:
                pass
    
        # Running Simlr 
        simlr = SIMLR.SIMLR_LARGE(num_of_rank=self.n_components, num_of_neighbor=self.n_neighbours, max_iter=self.max_iter)
        S, F,val, ind = simlr.fit(X)
        self.results = F
 
class Phenograph(UnsupervisedMethods):

    def __init__(self, n_neighbours, threads):
        UnsupervisedMethods.__init__(self)
        self.n_neighbours = n_neighbours
        self.threads = threads

    def apply(self):

        communities, graph, Q = phenograph.cluster(data=self.matrix,k=self.n_neighbours,n_jobs=threads)
        self.results = communities




















        

