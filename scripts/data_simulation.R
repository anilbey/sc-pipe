suppressMessages(library(splatter))
suppressMessages(library(rhdf5))

splat_simulate = function(data, facLoc=1, facScale=0.3, deProb=0.1, dropoutPresent=FALSE, groupProb=c(0.33,0.33,0.34))
{
   # Estimate parameters from the data (it takes time, could be moved out of the function)
    params = splatEstimate(data)
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
    h5write(gene_names, f_name, "gene_names")
    h5write(simulated_matrix, f_name,"matrix")
    h5write(cell_names, f_name, "cell_names")
    h5write(cell_groups, f_name, "cell_groups")
    h5write(FALSE, f_name, "cells_on_rows")
}

data_simulation = function (input_file, output_file, group_prob, dropout_present,  loc_factor, threads)
{
    # R uses column-major order therefore the HDF5 matrix is transposed
    # automatically. Python and hdf5 use row-major order
    print('dropout present value:')
    print(dropout_present)
    # loc_factor is a string since we read it from the wildcards!
    loc_factor = as.numeric(loc_factor)

    h5f = H5Fopen(input_file)
    data = h5f$matrix
    H5Fclose(h5f)
    print('*******shape of the matrix****')
    print(dim(data))
    sg = splat_simulate(data = data, groupProb = group_prob, dropoutPresent=dropout_present,  facLoc = loc_factor)
    write_hdf5(sg, output_file)

}

# snakemake@input is a list
data_simulation(snakemake@input[[1]], snakemake@output[[1]], snakemake@params[['group_prob']], snakemake@params[['dropout_present']],  snakemake@wildcards[['loc']], snakemake@threads)
