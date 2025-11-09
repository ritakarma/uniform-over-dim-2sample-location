#load the required files from the same directory
source("functions.R")

# Load required packages
library("ggplot2")
library("gridExtra")
# Colon tissue data is obtained from the package "plsgenomics".
# The dataset is also available at -- http://genomics-pubs.princeton.edu/oncology/affydata/index.html
library("plsgenomics")


# Loading Colon tissue data
data(Colon)
dataX = Colon$X[Colon$Y == 2, ]
dataY = Colon$X[Colon$Y == 1, ]
rm(Colon)

# Function to perform all high dimensional tests
all_tests <- function(dataX, dataY, vec_beta = c(0.25), nsim = 10000) {
  ##----------------------------------------------------------------------------
  ## Input: 
  ##   dataX: n1 * p matrix -rows are X observations  
  ##   dataY: n2 * p matrix -rows are Y observations
  ##   nsim: number of simulations for cut-off calculation
  ##   vec_beta:  vector specifying beta parameter values for second estimator
  ## Output: 
  ##   a vector of p-values of different tests
  ##----------------------------------------------------------------------------
  
  ## Performing our tests
  vec1 = perform_test(dataX = dataX, dataY = dataY, h = h_diff, nsim = 10000, vec_beta = vec_beta)
  vec2 = perform_test(dataX = dataX, dataY = dataY, h = h_spatial, nsim = 10000, vec_beta = vec_beta)
  vec = c(vec1, vec2)  
  
  if (length(beta) > 0) {
    names(vec) = c("difplain", paste0("diftap_", vec_beta), "spaplain", paste0("spatap_", vec_beta))
  }
  else {
    names(vec) = c("difplain", "spaplain")
  }
  
  ## Performing other tests
  pval = HDNRA::ZGZC2020.TS.2cNRT(dataX, dataY)$p.value # Zhang et al 2020 Jasa (ZGZC2020)
  vec = c(vec, ZGZC2020 = as.numeric(pval))   
  
  pval = highmean::apval_Bai1996(dataX, dataY)$pval # Bai and Saranadasa 1996 (BS1996)
  vec = c(vec, BS1996 = as.numeric(pval))
  
  pval = highmean::apval_Cai2014(dataX, dataY)$pval # Cai et al 2014 (CLX2014)
  vec = c(vec, CLX2014 = as.numeric(pval))
  
  pval = highmean::apval_Chen2010(dataX, dataY)$pval # Chen and Qin 2010 (CQ2010)
  vec = c(vec, CQ2010 = as.numeric(pval))
  
  pval = highmean::apval_Chen2014(dataX, dataY)$pval # Chen et al 2014 (CLZ2014)
  vec = c(vec, CLZ2014 = as.numeric(pval))
  
  pval = highmean::apval_Sri2008(dataX, dataY)$pval # Srivastava and Du 2008 (SD2008)
  vec = c(vec, SD2008 = as.numeric(pval))
  
  return(vec)
  
}


result = all_tests(dataX = dataX, dataY = dataY, nsim = 100000)
result_block = c()
for (k in 1:50) {
  result_block = rbind(result_block, all_tests(dataX[ , (40 * (k-1) + 1):(40 * k)], dataY[ , (40 * (k-1) + 1):(40 * k)]))
}
result_block_avg = colSums(result_block) / nrow(result_block)

print("p-value of different tests for the whole dataset")
print(result)

print("average p-value of different tests corresponding to the data blocks")
print(result_block_avg)

print("The tests 'difplain', 'diftap_0.25', 'spaplain' and 'spatap_0.25' correspond to 'KCDG2025^1', 'KCDG2025^2', 'sKCDG2025^1', and 'sKCDG2025^2' respectively.")


# Prepare test names and indices for histogram
test_names <- colnames(result_block)

# Changing test names: "difplain" - "KCDG2025^1", "diftap_0.25" - "KCDG2025^2", "spaplain" - "sKCDG2025^1", "spatap_0.25" - "sKCDG2025^2" 

test_names[1:4] <- c("KCDG2025^1", "KCDG2025^2", "sKCDG2025^1", "sKCDG2025^2")
tests <- c(5, 6, 7, 8, 10, 1, 2, 3, 4)  # total 9 histograms

# Create a list of ggplot histograms
plots <- lapply(seq_along(tests), function(k) {
  i <- tests[k]
  ggplot(data.frame(p = result_block[, i]), aes(x = p)) +
    geom_histogram(
      breaks = seq(0, 1, by = 0.02),
      fill = "lightgray",
      color = "black"
    ) +
    labs(
      x = test_names[i],
      y = "Frequency"
    ) +
    theme_minimal(base_size = 26) +  # increase overall base text size
    theme(
      plot.margin = margin(10, 10, 10, 10),
      axis.title.x = element_text(size = 30, margin = margin(t = 12)),
      axis.title.y = element_text(size = 30, margin = margin(r = 12)),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 22),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
})

# Arrange in 2 rows: 5 on first row, 4 centered on second
row1 <- do.call(gridExtra::arrangeGrob, c(plots[1:5], nrow = 1))
row2 <- do.call(gridExtra::arrangeGrob, c(plots[6:9], nrow = 1))

# Add blank space on both sides of the second row to center it
blank <- grid::rectGrob(gp = grid::gpar(col = NA))
row2_centered <- gridExtra::arrangeGrob(blank, row2, blank, nrow = 1, widths = c(1, 8, 1))

# Combine both rows vertically
final_plot <- gridExtra::arrangeGrob(row1, row2_centered, ncol = 1, heights = c(1, 1))


# Save as PDF
ggsave("colon_data_histogram.pdf", final_plot, width = 32, height = 16)



  