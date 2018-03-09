suppressMessages(library(splatter))
suppressMessages(library(rhdf5))
library(argparse)

parser <- ArgumentParser(description= 'simulates data')
parser$add_argument('--input', help= 'rda input containing the estimated parameters' , required=TRUE)
parser$add_argument('--group_prob', type='double', nargs='*', required=TRUE)
parser$add_argument('--dropout_present', type='logical', required=TRUE)
parser$add_argument('--output', required=TRUE)
parser$add_argument('--loc', required=TRUE, type='double')
parser$add_argument('--de_prob', required=TRUE, type='double')
args = parser$parse_args()

print('dropout present value:')
print(args$dropout_present)

splat_simulate = function(input_rda, facLoc=1, facScale=0.3, deProb=0.1, dropoutPresent=FALSE, groupProb)
{
   # Estimate parameters from the data (it takes time, could be moved out of the function)
    # loads the params object from the rda
    load(input_rda)
    params = setParam(params, "group.prob", groupProb)
    params = setParam(params, "de.facScale", facScale)
    params = setParam(params, "de.prob", deProb)
    params = setParam(params, "de.facLoc", facLoc)
    params = setParam(params, "dropout.present", dropoutPresent)
    print(params)
    sg = splatSimulateGroups(param=params, verbose=TRUE)
    return (sg)
}

write_hdf5 = function(sg, f_name)
{
    print('write hdf5 function')
    print(f_name)
    simulated_matrix = assays(sg)$counts
    gene_names = rownames(simulated_matrix)
    cell_names = colnames(simulated_matrix)
    cell_groups = sg$Group
    h5createFile(f_name)
    h5createGroup(f_name, "cell_attrs")
    h5createGroup(f_name, "gene_attrs")
    h5write(gene_names, f_name, "gene_attrs/gene_names")
    h5write(simulated_matrix, f_name,"matrix")
    h5write(cell_names, f_name, "cell_attrs/cell_names")
    h5write(cell_groups, f_name, "cell_attrs/cell_groups")
    h5write(FALSE, f_name, "cell_attrs/cells_on_rows")

    X = colnames(rowData(sg))
    for (deg in X[grepl("^DEFacGroup*", X)])
    {
        deg_values = rowData(sg)[deg][[1]]
        h5write(deg_values, f_name, paste("gene_attrs/",deg ,sep = ''))
    }
    # closes all open HDF5 handles in the environment 
    H5close()

}

data_simulation = function (input_rda, output_file, group_prob, dropout_present,  loc_factor, de_prob)
{
    # R uses column-major order therefore the HDF5 matrix is transposed
    # automatically. Python and hdf5 use row-major order

    # loc_factor is a string since we read it from the wildcards!
    loc_factor = as.numeric(loc_factor)

    sg = splat_simulate(input_rda = input_rda, groupProb = group_prob, dropoutPresent=dropout_present,  facLoc = loc_factor, deProb=de_prob)
    write_hdf5(sg, output_file)

}

data_simulation(args$input, args$output, args$group_prob, args$dropout_present,  args$loc, args$de_prob)
