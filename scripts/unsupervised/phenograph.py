import phenograph
from .unsupervised_method import UnsupervisedMethod

class Phenograph(UnsupervisedMethod):

    def __init__(self, n_neighbours, threads):
        super().__init__()
        self.n_neighbours = n_neighbours
        self.threads = threads

    def apply(self):
        communities, graph, Q = phenograph.cluster(data=self.matrix,k=self.n_neighbours,n_jobs=self.threads)
        self.results = communities
