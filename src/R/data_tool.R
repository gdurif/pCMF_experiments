## functions to manage single-cell RNAseq data

## Delete genes that are expressed (count >= 'min_count') in more than 'n_cell'
## If 'percentage' is TRUE, then 'n_cell' is a proportion
## The check can be done depending on a partition of cells (according to 'condition' vector) 
del_zero_ind <- function(data, n_cell = 0, 
                         min_count = 0, 
                         percentage = FALSE, 
                         conditions) {
    x_col <- rep(FALSE, ncol(data))
    if(missing(conditions)){
        if(percentage){
            x_col <- colSums(data > min_count) / nrow(data) > n_cell & colSums(data > min_count) > 1
        } else {
            x_col <- colSums(data > min_count) > n_cell
        }
    } else {
        conditions <- as.data.frame(conditions)
        for(i in ncol(conditions)){
            for(state in levels(as.factor(as.vector(conditions[,i])))){
                r_select <- conditions[,i] == state
                if(percentage){
                    x_col_state <- colSums(data[r_select,] > min_count) / nrow(data[r_select,]) > n_cell & colSums(data > min_count) > 1
                } else {
                    x_col_state <- colSums(data[r_select,] > min_count) > n_cell
                }
                x_col <- ifelse(x_col_state, TRUE, x_col)
            }
        }
    }
    return(x_col)
}