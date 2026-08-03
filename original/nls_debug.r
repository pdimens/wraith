#! /usr/bin/env R

simulateWrath <- function(data){
    # assign NAs to anything below and including the diagonal
    #data_m <- as.matrix(data)
    data[lower.tri(data, diag = TRUE)] <- NA

    # matrix distribution taking into account position
    data_df <- as.data.frame(data)

    # get row numbers and save them as a column
    data_df$nrow <- seq.int(nrow(data_df))

    # get column numbers and save them as a column
    data_df <- pivot_longer(data_df, !nrow, names_to = "ncol")

    # omit NAs
    data_df <- na.omit(data_df)

    # modify ncols to make them integers
    data_df$ncol <- as.integer(gsub("V", "", data_df$ncol))

    # calculate the index by the distance to the matrix
    data_df_upper <- group_by(data_df, nrow) %>%
        mutate(index = ncol - nrow)

    # create dataframe with x and y values to calculate z-scores from
    points <- data_df_upper %>%
        ungroup() %>%
        dplyr::select(value, index, ncol, nrow)
    colnames(points) <- c("y", "x", "ncol", "nrow")

    #### PART 2: Z-SCORE CALCULATION ####
    # calculate z scores grouping values by their distance to the diagonal
    scaled_full_df <- tibble(z_score = numeric(), y = numeric(), x = numeric(), nrow = numeric(), ncol = numeric())

    for (index in unique(points$x)) {
        group_indices <- which(points$x == index)
        group_values <- points[group_indices, ]$y
        scaled_group <- scale(na.omit(group_values))
        scaled_df <- tibble(z_score = scaled_group[, 1], y = points[group_indices, ]$y, x = index, nrow = points[group_indices, ]$nrow, ncol = points[group_indices, ]$ncol)
        scaled_full_df <- full_join(scaled_full_df, scaled_df, by = c("z_score", "y", "x", "nrow", "ncol"))
    }
    return(list(fulldf = scaled_full_df, points = points))
}

### NOTES:
## the points dataframe, which looks like:by$points
# A tibble: 4,950 × 4
#       y     x  ncol  nrow
#   <dbl> <int> <int> <int>
# 1 11.4      1     2     1
# 2 12.9      2     3     1
# 3 16.7      3     4     1
# 4 17.4      4     5     1
# 5 16.8      5     6     1
# 6 13.3      6     7     1
# 7 16.9      7     8     1
# 8  6.30     8     9     1
# 9  1.12     9    10     1
#10  7.29    10    11     1

## COLUMNS:
# y = jaccard score
# x = index along diagonal
# ncol = column index of y (and Zscore in scaled_full_df)
# nrow = like ncol but row index