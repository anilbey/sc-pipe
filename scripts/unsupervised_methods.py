from abc import ABCMeta, abstractmethod
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.decomposition import FactorAnalysis as sk_fa
import h5py
from ZIFA import ZIFA, block_ZIFA
import numpy as np
import pandas as pd
import SIMLR
import phenograph
from sklearn.metrics.cluster import silhouette_score

class UnsupervisedMethod(metaclass=ABCMeta):

    def __init__(self):
        self.matrix = None
        self.barcodes = None
        self.results = None

    def load_from_csv(self, input_file):
        df = pd.read_csv(input_file, header=None)
        self.barcodes = df[0].values
        df = df.drop(axis=1, columns=[0])
        self.matrix = df.values

    def load_from_hdf5(self, input_file):
        h5f = h5py.File(input_file, 'r')
        self.matrix = h5f['matrix'].value
        barcodes = h5f['cell_attrs']['cell_names'].value
        decoder = np.vectorize(lambda t: t.decode('UTF-8'))
        self.barcodes = decoder(barcodes)
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

class Pca(UnsupervisedMethod):

    def __init__(self, n_components):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components

    def apply(self):
        pca = PCA(n_components = self.n_components)
        self.results = pca.fit_transform(self.matrix)

class Tsne(UnsupervisedMethod):

    def __init__(self, n_components, init):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components
        self.init = init

    def apply(self):
        tsne = TSNE(init=self.init, n_components=self.n_components)
        self.results = tsne.fit_transform(self.matrix)

class FactorAnalysis(UnsupervisedMethod):

    def __init__(self, n_components):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components

    def apply(self):
        fa = sk_fa(n_components = self.n_components)
        self.results = fa.fit_transform(self.matrix)

class Zifa(UnsupervisedMethod):

    def __init__(self, n_components, n_blocks):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components
        self.n_blocks = n_blocks

    def apply(self):
        Zhat, params = block_ZIFA.fitModel(self.matrix, self.n_components, n_blocks = self.n_blocks)
        self.results = Zhat

class Simlr(UnsupervisedMethod):
    
    def __init__(self,n_components, pca_components, n_neighbours, max_iter, threads):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components
        self.pca_components = pca_components
        self.n_neighbours = n_neighbours
        self.max_iter = max_iter
        self.threads = threads

    def set_params(self, n_clusters):
        self.n_components = n_clusters
        return self

    def fit_predict(self, matrix):
        X = matrix
    
        # Selecting 500 features with PCA
        if X.shape[1]>500:
        # fast_pca assumes the number of cells > 500 therefore try-catch
            try:
                X = SIMLR.helper.fast_pca(X,self.pca_components)
            except:
                pass    
        # Running Simlr 
        simlr = SIMLR.SIMLR_LARGE(num_of_rank=self.n_components, num_of_neighbor=self.n_neighbours, max_iter=self.max_iter)
        S, F, val, ind = simlr.fit(X)
        return F

    def apply(self):
        self.results = fit_predict(self.matrix)

class Phenograph(UnsupervisedMethod):

    def __init__(self, n_neighbours, threads):
        super().__init__()
        self.n_neighbours = n_neighbours
        self.threads = threads

    def apply(self):
        communities, graph, Q = phenograph.cluster(data=self.matrix,k=self.n_neighbours,n_jobs=self.threads)
        self.results = communities

class Silhouette(UnsupervisedMethod):

    def __init__(self, model, k_min, k_max, metric):
        UnsupervisedMethod.__init__(self)
        self.model = model
        self.k_min = k_min
        self.k_max = k_max
        self.metric = metric
        

    @classmethod
    def init_with_kmeans(cls, n_init, n_jobs, k_min, k_max, metric):
       
        # K-means falls in local minima. That’s why it can be useful to restart it several times using n_init
        # The classical EM-style algorithm is "full" and it is suggested for sparse data
        
        model = KMeans(algorithm='full', n_init=n_init, n_jobs=n_jobs)
        return cls(model, k_min, k_max, metric)

    @classmethod
    def init_with_hierarchical(cls, h_affinity, h_linkage, k_min, k_max, metric):
        model = AgglomerativeClustering(affinity=h_affinity, linkage=h_linkage, compute_full_tree=False)
        return cls(model, k_min, k_max, metric)
    
    @classmethod
    def init_with_simlr(cls, n_components, pca_components, n_neighbours, max_iter, threads, k_min, k_max, metric):
        model = Simlr(n_components, pca_components, n_neighbours, max_iter, threads)
        return cls(model, k_min, k_max, metric)

    def apply(self):
        
        k_range = range(self.k_min, self.k_max)
        if isinstance(self.model, Simlr):
            predicted_labels =  [KMeans().fit_predict(self.model.set_params(n_clusters=k).fit_predict(self.matrix)) for k in k_range]
        else:
            predicted_labels = [self.model.set_params(n_clusters=k).fit_predict(self.matrix) for k in k_range]

        silhouette_scores = [silhouette_score(X=self.matrix, labels=obj, metric=self.metric) for obj in predicted_labels]
        max_index = np.argmax(silhouette_scores)
        self.results = predicted_labels[max_index]









