##################################################
## File name: select_protein_coding_genes.R
## Author: Dr. Franziska Singer and Mustafa Anil Tuncel
## Date created: 22.03.2018
## R Version: 3.4.1
##################################################

library(argparse)
library(rhdf5)
library(biomaRt)


parser =  ArgumentParser(description= 'filters out non-coding genes')
parser$add_argument('--input', help= 'the input hdf5 file' , required=TRUE)
parser$add_argument('--output', help= 'the output hdf5 file containing the filtered results', required=TRUE)
args = parser$parse_args()

path = args$input
output_file = args$output

select_protein_coding_genes = function(path, output_file)
{
  h5f = H5Fopen(path)

  # read the gene ids
  gene_ids = h5f$gene_attrs$gene_ids

  # map ensembl to entrez gene ids, only take into account protein coding genes
  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

  entrezGeneMapping_proteinCoding = getBM(attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"),
      filters = c("ensembl_gene_id","biotype"), values= list(gene_ids,'protein_coding'), mart= mart, uniqueRows=T)

  # filter NAs and multiple entries
  entrezGeneMapping_proteinCoding_noNa = na.omit(entrezGeneMapping_proteinCoding)
  entrezGeneMapping_proteinCoding_noNa_unique = entrezGeneMapping_proteinCoding_noNa[!duplicated(entrezGeneMapping_proteinCoding_noNa$ensembl_gene_id),]

  print(paste('Mapped to entrez and kept only protein coding genes. Remaining #genes:', dim(entrezGeneMapping_proteinCoding_noNa_unique)[1], sep =' '))

  length(entrezGeneMapping_proteinCoding_noNa_unique$ensembl_gene_id)

  # creates a boolean (logical) mask vector
  res = gene_ids %in% entrezGeneMapping_proteinCoding_noNa_unique$ensembl_gene_id

  # retrieve the rest of the hdf5 attributes and apply the mask
  data = h5f$matrix
  f_name = 'coding_genes_output_hdf5'
  gene_names = h5f$gene_attrs$gene_names[res]
  gene_ids = h5f$gene_attrs$gene_ids[res]
  cell_names = h5f$cell_attrs$cell_names
  filtered_matrix = data[res,]

  H5Fclose(h5f)

  f_name = output_file

  # write the results
  h5createFile(f_name)
  h5createGroup(f_name, "cell_attrs")
  h5createGroup(f_name, "gene_attrs")
  h5write(gene_names, f_name, "gene_attrs/gene_names")
  h5write(gene_ids, f_name, "gene_attrs/gene_ids")
  h5write(filtered_matrix, f_name,"matrix")
  h5write(cell_names, f_name, "cell_attrs/cell_names")
}

select_protein_coding_genes(path,output_file)
