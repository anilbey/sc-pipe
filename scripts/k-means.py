from sklearn.cluster import KMeans
import numpy as np



def apply_kmeans(input_file, output_file, n_clusters, threads):
    
    data = np.genfromtxt(input_file, delimiter=',')
    y_pred = KMeans(n_clusters=n_clusters, random_state=123).fit_predict(data)
    np.savetxt(output_file, y_pred, fmt='%i', delimiter=",")

apply_kmeans(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_clusters, snakemake.threads)