library(dplyr)
library(ggplot2)
library(tidyr)

# define function for impulse response visualisation
irf_plot_replication <- function(data){
  
  # Compute quantiles for each response variable and quarter
  quantiles <- apply(data,c(2,3),quantile,probs=c(0.01,0.25,0.5,0.75,0.99))
  
  # Transform quantile array into table for visualising
  viz_df <- as.data.frame.table(quantiles) %>% 
    # Proper renaming
    rename(c("measure"="Var1","indic"="Var2","quarter"="Var3","value"="Freq")) %>% 
    mutate(measure = as.character(measure),
           indic = as.character(indic), # Country-Variable response indicator
           quarter = as.numeric(quarter), # time axisi numerical from factor
           value = value * 100 # Plots shows percent response, so multiply by 100
           ) %>% 
    # Separate country-variable indicator into components
    separate(indic, into = c("country", "response_var"), sep = "\\.") %>% 
    # Create grouping variable on regions as indicated in paper
    mutate(region = case_when(
      country %in% c("EA","UK","CA","AU","NZ","CH","NO","SE","DK","IS") ~ "Advanced Economies",
      country %in% c("CZ","HU","PL","SK","SI","BG","RO","HR","AL","RS","TR",
                     "LT","LV","EE","RU","UA","BY","GE") ~ "Emerging Europe",
      country %in% c("CN","KR","JP","PH","SG","TH","ID","IN","MY") ~ "Asia",
      country %in% c("AR","BR","CL","MX","PE") ~ "Latin America",
      TRUE ~ NA_character_),
      # Rename quantile names in column format
      measure = case_when(
        measure == "1%" ~ "p1",
        measure == "25%" ~ "p25",
        measure == "50%" ~ "p50",
        measure == "75%" ~ "p75",
        measure == "99%" ~ "p99"
      )) %>% 
    select(-country) %>% 
    # Deselect US country and long-term rates, trade balance variables
    filter(!is.na(region) & response_var %in% c("y","Dp","stir", "rer")) %>% 
    # Group by relavant variables and compute regional mean response
    group_by(measure, region, response_var, quarter) %>% 
    summarise(value = mean(value,na.rm=TRUE)) %>% 
    ungroup() %>% 
    # Make table wide so each quantile is a column
    pivot_wider(names_from=measure, values_from = value)
  
  # Define variable labels
  var_labels <- c("y" = "Real GDP","Dp" = "Inflation",
                  "stir" = "Short-term rates", "rer"="Real effective ER")
  
  # Start making plot
  plot <- viz_df %>% 
    # Assign factor structure to mirror plot layout in paper
    mutate(response_var = factor(response_var, levels = c("y","Dp","stir", "rer"))) %>% 
    # Plot median, P25 and P75 response 
    ggplot(aes(x=quarter, y=p50, ymin=p25, ymax=p75)) +
    geom_hline(yintercept = 0, color="red") +
    geom_ribbon(fill="grey", alpha=0.2) +
    geom_line() +
    # grid of plots for region rows and variable columns
    facet_grid(region ~ response_var,
               labeller = labeller(response_var = var_labels),
               scales = "free_y") + 
    # Rest for visual beautification
    theme_light() +
    ylab("")+
    xlab("Quarters") +
    theme(plot.title = element_text(size = 11, hjust=0.5),
          axis.title.y = element_text(size=11))
  
  return(plot)
}

#IRF_post is the posterior of structural IRF, 1st dim: Draws, 2nd dim: Response, 3rd.dim: shock, 4rd. dim: horizon
irf_data <- readRDS("data/BGVAR_IRF_estimates.rds")
dimnames(irf_data)

# Make plots for the three reported shocks (poil is disregarded)
yshock_plot <- irf_plot_replication(irf_data[,,1,])

pshock_plot <- irf_plot_replication(irf_data[,,2,])

rshock_plot <- irf_plot_replication(irf_data[,,3,])

# Save plots 
ggsave(plot= yshock_plot, filename = "data/yshock_plot.png")
ggsave(plot= pshock_plot, filename = "data/pshock_plot.png")
ggsave(plot= rshock_plot, filename = "data/rshock_plot.png")