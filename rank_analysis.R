#FINAL ANALYSIS FAMILY
print("welcome")

#-------LIBRARIES-------
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse) #To parse the biome names
library(forcats)  # Needed for fct_reorder()
library(gridExtra) #to plot multiple graphs
library(stringr) #to modify categories
library(rlang) #use as_name()
#library(kableExtra) #for correlation tables
#library(indicspecies) # For species specificity
#library(data.table) #to calculate specificity of bgc

#-------PREPROCESSING FUNCTIONS-------

#Function to merge regions taxonomy and proper name from taxdump
name_regions <- function(regions, taxdump, rank_id_col = "genus_id", rank_name_col = "genus_name") {
  # Merge regions_taxonomy with gen_taxdump to get genus_id
  regions_with_taxon <- regions %>%
    inner_join(taxdump, by = "tax_id")  # Join by tax_id
  
  # Now, extract the name of each genus_id by finding where genus_id matches tax_id in gen_taxdump
  taxon_names <- taxdump %>%
    select(tax_id, name) %>%  # Extract only relevant columns
    rename(!!rank_id_col := tax_id, !!rank_name_col := name)  # Rename dynamically
  
  # Merge the genus name into the dataset
  regions_with_taxon_named <- regions_with_taxon %>%
    left_join(taxon_names, by = rank_id_col)  # Match genus_id with its name
  
  # Replace empty string with NA
  regions_with_taxon_named[regions_with_taxon_named == ''] <- NA
  
  return(regions_with_taxon_named)
}

#Changes labels that account for less than the treshold to the most common one, considered as mislabelings.
adjust_labels <- function(data, id_col, category_col, threshold = 5) {
  # Extract unique pairs of id and category
  id_to_category <- data %>%
    select({{id_col}}, {{category_col}}) %>%
    distinct()
  
  # Vector of IDs with more than one category label
  ids_with_multiple_labels <- names(which(table(pull(id_to_category, {{id_col}})) > 1))
  
  # Create a duplicate column to adjust the labels
  adjusted_col_name <- paste0("adjusted_", as_label(ensym(category_col)))
  data <- data %>%
    mutate(!!adjusted_col_name := {{category_col}})
  
  # Loop through each ID and adjust the labels
  for (id in ids_with_multiple_labels) {
    id <- as.character(id)
    
    # Filter data for the current ID
    filtered_data <- data %>%
      filter({{id_col}} == id)
    
    # Count occurrences of each category
    category_counts <- filtered_data %>%
      count({{category_col}}, name = "count") %>%
      mutate(percentage = count / sum(count) * 100)
    
    # Identify labels below the threshold
    rare_labels <- category_counts %>%
      filter(percentage < threshold) %>%
      pull({{category_col}})
    
    # Skip to next ID if no labels are below the threshold
    if (length(rare_labels) == 0) next
    
    # Find the most common category
    most_common_label <- category_counts %>%
      arrange(desc(count)) %>%
      slice(1) %>%
      pull({{category_col}})
    
    # Adjust the labels in the main dataframe
    data <- data %>%
      mutate(!!adjusted_col_name := ifelse(
        {{id_col}} == id & {{category_col}} %in% rare_labels,
        most_common_label,
        .data[[adjusted_col_name]]
      ))
  }
  
  return(data)
}

#Splits longest_biome into its parts
adjust_biome <- function(data){
  
  #Replace NA in longest_biome to root
  data[["longest_biome"]][is.na(data[["longest_biome"]])] <- "root"
  # Split the column and create the new columns
  data <- data %>%
    mutate(
      split_values = str_split(longest_biome, ":", simplify = FALSE),
      superbiome = map_chr(split_values, ~ .x[2] %||% NA_character_),
      biome = map_chr(split_values, ~ .x[3] %||% NA_char),
      subbiome = map_chr(split_values, ~ if (length(.x) >= 4) paste(.x[3:4], collapse = ":") else NA),
      infrabiome = map_chr(split_values, ~ if (length(.x) >= 5) paste(.x[3:5], collapse = ":") else NA),
      microbiome = map_chr(split_values, ~ if (length(.x) >= 6) paste(.x[3:6], collapse = ":") else NA)
    ) %>% 
    select(-split_values) # Drop the intermediate split_values column
  
  return(data)
}

#-------DATA PREPROCESSING-----------

#Read taxonomy regions file
regions_taxonomy <- read.delim(file = "atlas_tables/regions_taxonomy.tsv", sep = "\t")
#regions_taxonomy <- na.omit(regions_taxonomy) #Gets filtered out later anyways

#Read family taxdump file
fam_taxdump<-read.delim(file = "atlas_tables/family_taxdump.tsv", sep = "\t", header=FALSE)
colnames(fam_taxdump) <- c("tax_id","parent","rank","name", "family_id")

#Read species taxdump file
spec_taxdump<-read.delim(file = "atlas_tables/species_taxdump.tsv", sep = "\t", header=FALSE)
colnames(spec_taxdump) <- c("tax_id","parent","rank","name", "species_id")

#Merging all taxdump files
#complete_taxdump <- full_join(fam_taxdump, gen_taxdump[c("tax_id","genus_id")], by = "tax_id")
#complete_taxdump <- full_join(complete_taxdump, spec_taxdump[c("tax_id","species_id")], by = "tax_id")      
#Write table
#write.table(complete_taxdump, file='complete_taxdump.tsv', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

#MERGING DATASETS
#Creating regions table with genus taxonomy
regions_with_fam_named <-  name_regions(regions_taxonomy, fam_taxdump, rank_id_col = "family_id", rank_name_col = "family_name")
#Creating regions table with genus taxonomy
regions_with_spec_named <-  name_regions(regions_taxonomy, spec_taxdump, rank_id_col = "species_id", rank_name_col = "species_name")

#PROCESSING LONGEST_BIOME
regions_with_fam_named <- adjust_biome(regions_with_fam_named)
#Filter out Mixed NA, and Engineered biomes
regions_with_fam_named <- regions_with_fam_named[regions_with_fam_named$superbiome %in% c("Host-associated", "Environmental"), ]
#Before: 1,670,503 After: 1,490,018, Lost:  10.8%
regions_with_spec_named  <- adjust_biome(regions_with_spec_named )
#Filter out Mixed NA, and Engineered biomes
regions_with_spec_named  <- regions_with_spec_named [regions_with_spec_named$superbiome %in% c("Host-associated", "Environmental"), ]
#Before: 777,406 , After: 703,811 Lost:  9.5%
#PROCESSING GCF IDs
#Turning gcf ids from integer to character
regions_with_fam_named[,14] <- sapply(regions_with_fam_named[,14], as.character)
#Turning gcf ids from integer to character
regions_with_spec_named[,14] <- sapply(regions_with_spec_named[,14], as.character)


#-------COUNT FUNCTIONS--------------- 

#Function for calculating region counts of chosen variable and representing it split by another variable
plot_bgc_split_counts <- function(variable, split_value, regions_with_genus_named, top = 0, plot = TRUE, split=TRUE, color = "orange2") {
  # Convert to symbols for tidy evaluation
  var_sym <- ensym(variable)
  split_sym <- ensym(split_value)
  
  # Count regions per variable
  variable_counts <- regions_with_genus_named %>%
    group_by(!!split_sym, !!var_sym) %>%  # Group by split_value and variable
    summarise(region_count = n(), .groups = "drop")
  
  # Subsetting variables with most regions
  top_variable_counts <- variable_counts %>% 
    filter(region_count >= top) %>%
    pull(!!var_sym)
  
  # Make counts for every split_value
  split_counts <- regions_with_genus_named %>%
    group_by(!!var_sym, !!split_sym) %>%
    tally() %>%
    pivot_wider(names_from = !!split_sym, values_from = n)
  
  
  # Replacing NA values with zero or "NA"
  split_counts <- split_counts %>%
    mutate(across(everything(), ~ if (is.numeric(.)) {
      replace(., is.na(.), 0)
    } else {
      replace(., is.na(.), "na")
    }))
  
  if(plot){
    # Selecting and plotting only variables with most regions
    top_split_counts <- split_counts %>%
      filter(!!var_sym %in% top_variable_counts)
    
    if(split){
      # Plotting
      plot <- top_split_counts %>%
        pivot_longer(-!!var_sym) %>%
        ggplot(aes(x = reorder(!!var_sym, -value), y = value, fill = name)) +
        geom_col(position = position_stack()) + 
        ylab("Number of Regions") + 
        xlab(paste(as.character(var_sym), "split by", as.character(split_sym) )) +
        theme(legend.position = "none", 
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }else{
      # Plotting
      plot <- top_split_counts %>%
        pivot_longer(-!!var_sym) %>%
        ggplot(aes(x = reorder(!!var_sym, -value), y = value)) +
        geom_col(position = position_stack(), fill = color, alpha = .8) + 
        ylab("Number of Regions") + 
        xlab(paste(as.character(var_sym) )) +
        theme(legend.position = "none", 
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    # Print the plot
    print(plot)
  }
  
  # Return the split counts
  return(split_counts)
}
#


#---------COUNT BIOME SPLIT BY FAMILY-----------------

#Counts of biome split by family
biome_family_counts <- plot_bgc_split_counts(biome, family_name,regions_with_fam_named, 0 , plot = FALSE)
subbiome_family_counts <- plot_bgc_split_counts(subbiome, family_name,regions_with_fam_named, 0, plot = FALSE)
infrabiome_family_counts <- plot_bgc_split_counts(infrabiome, family_name,regions_with_fam_named, 0, plot = FALSE)
microbiome_family_counts <- plot_bgc_split_counts(microbiome, family_name,regions_with_fam_named, 0, plot = FALSE)
#

#---------COUNT BIOME SPLIT BY SPECIES-----------------

#Counts of biome split by family
biome_spec_counts <- plot_bgc_split_counts(biome, species_name,regions_with_spec_named, 0 , plot = FALSE)
subbiome_spec_counts <- plot_bgc_split_counts(subbiome, species_name,regions_with_spec_named, 0, plot = FALSE)
infrabiome_spec_counts <- plot_bgc_split_counts(infrabiome, species_name,regions_with_spec_named, 0, plot = FALSE)
microbiome_spec_counts <- plot_bgc_split_counts(microbiome, species_name,regions_with_spec_named, 0, plot = FALSE)
#

#-------DIVERSITY FUNCTIONS-----------

#Function that takes a column  and calculates its simpson index
simpson_idx <- function(column){
  N = sum(column)
  sum = 0
  for (row in column){
    n = row
    sum = sum + ((n * (n - 1))/(N * (N - 1 )))
  }
  if (is.nan(sum)){
    sum = 1
  }
  return(sum)
}

#Function that reverses the simpson index
reverse_simpson <- function(column){
  return(1 - simpson_idx(column))
}

#Function to calculate shannon-wiener index for column
shannon_wieder_idx <- function(column) {
  N = sum(column)
  sum = 0
  for (row in column){
    n = row
    p = n/N
    h = (p * log(p))
    if (is.nan(h)){
      h = 0
    }
    sum = sum + h
  }
  return(-sum)
}

#Calculate specificity values of a variable for a given split value
plot_split_diversity <- function(counts, variable, split_value ,min_richness = 0, plot = TRUE) {
  # Convert to symbol for tidy evaluation
  var_sym <- ensym(variable)
  
  # Calculating the number of different split_values found in each variable
  diversity <- data.frame(
    variable = counts[[as.character(var_sym)]], 
    richness = rowSums(counts[, -1] > 0)
  )
  #Changing variable column name
  colnames(diversity)[1] <- as.character(var_sym) 
  
  # Transposing table for easier calculations
  t_counts <- setNames(
    data.frame(t(counts[,-1])), 
    counts[[as.character(var_sym)]]
  )
  
  # Applying function to every column and adding to type_diversity table
  diversity$r_simpson <- sapply(t_counts, reverse_simpson)
  diversity$shannon <- sapply(t_counts, shannon_wieder_idx)
  
  if(plot){
    
    # Filter by minimum richness
    filtered_diversity <- diversity %>%
      filter(richness > min_richness)
    
    # Plotting
    plot <- ggplot(
      data = filtered_diversity %>%
        mutate(!!var_sym := fct_reorder(!!var_sym, shannon, .desc = TRUE)) %>%
        gather(Variable, value, -c(!!var_sym, richness)),  # Include r_simpson and shannon
      aes(x = !!var_sym, y = value, fill = Variable)
    ) + 
      geom_bar(stat = 'identity', position = 'dodge') +
      #geom_hline(yintercept = mean(filtered_diversity$shannon), color = "skyblue") +
      #geom_hline(yintercept = mean(filtered_diversity$r_simpson), color = "red3") +
      geom_hline(yintercept = 3, color = "skyblue") +
      geom_hline(yintercept = 0.65, color = "red3") +
      labs(y = paste("Diversity of ", split_value ), x = as.character(var_sym)) +
      scale_y_continuous(breaks=seq(0,5.6,0.2)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # Print the plot
    print(plot)
  }
  
  # Return type_diversity
  return(diversity)
}

#Calculating diversity individualy per biome
multiple_diversity <- function(data, column) {
  
  # Initialize an empty list to store resulting dataframes
  results <- list()
  # Loop through each unique value in the specified column
  
  for (value in unique(data[[column]])) {
    if(!is.na(value)) {
      print(value)
      # Filter the data for the current unique value
      subset_data <- data %>% filter(.data[[column]] == value)
      
      # Run the specified functions
      genus_name <- value  # Assuming genus_name is the unique value of the column
      gcf_genus_counts <- plot_bgc_split_counts(genus_name, bigslice_gcf_id, subset_data, plot=FALSE)
      final_df <- plot_split_diversity(gcf_genus_counts, genus_name, "gcf", plot =FALSE)
      
      # Add a column with the unique value of the specified column
      final_df[[column]] <- value
      
      # Store the dataframe in the results list
      results[[as.character(value)]] <- final_df
    }
  }
  
  # Combine all dataframes by row
  combined_results <- bind_rows(results)
  
  return(combined_results)
}


#---------FAMILY DIVERSITY FOR BIOMES---------

#Diversity values of genus in biome
biome_fam_diversity <- plot_split_diversity(biome_family_counts, biome, "family", plot = FALSE)
subbiome_fam_diversity <- plot_split_diversity(subbiome_family_counts, subbiome, "family", plot = FALSE)
infrabiome_fam_diversity <- plot_split_diversity(infrabiome_family_counts, infrabiome, "family", plot = FALSE)
microbiome_fam_diversity <- plot_split_diversity(microbiome_family_counts, microbiome, "family",plot = FALSE)

#---------SPECIES DIVERSITY FOR BIOMES---------

#Diversity values of genus in biome
biome_spec_diversity <- plot_split_diversity(biome_spec_counts, biome, "species", plot = FALSE)
subbiome_spec_diversity <- plot_split_diversity(subbiome_spec_counts, subbiome, "species", plot = FALSE)
infrabiome_spec_diversity <- plot_split_diversity(infrabiome_spec_counts, infrabiome, "species", plot = FALSE)
microbiome_spec_diversity <- plot_split_diversity(microbiome_spec_counts, microbiome, "species",plot = FALSE)

#---------COMPARING FAMILY DIVERSITY FOR BIOME LEVELS--------

#COMPARING SHANNON
# Create a combined dataframe
biome_level_family_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_fam_diversity$shannon),
  data.frame(Level = "Subbiome", Shannon = subbiome_fam_diversity$shannon),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_fam_diversity$shannon),
  data.frame(Level = "Microbiome", Shannon = microbiome_fam_diversity$shannon)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

# Create the boxplot
ggplot(biome_level_family_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Genus Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_fam_diversity$shannon,microbiome_fam_diversity$shannon)
# Kolmogorov-Smirnov Test: D = 0.25455, p-value = 0.4388
wilcox.test(biome_fam_diversity$shannon,microbiome_fam_diversity$shannon)
# Wilcoxon Rank-Sum Test: W = 216, p-value = 0.4947
# Test indicate distributions are not significantly different

#---------COMPARING SPECIES DIVERSITY FOR BIOME LEVELS--------

#COMPARING SHANNON
# Create a combined dataframe
biome_level_species_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_spec_diversity$shannon),
  data.frame(Level = "Subbiome", Shannon = subbiome_spec_diversity$shannon),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_spec_diversity$shannon),
  data.frame(Level = "Microbiome", Shannon = microbiome_spec_diversity$shannon)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

# Create the boxplot
ggplot(biome_level_species_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Genus Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_spec_diversity$shannon,microbiome_spec_diversity$shannon)
# Kolmogorov-Smirnov Test: 
wilcox.test(biome_spec_diversity$shannon,microbiome_spec_diversity$shannon)
# Wilcoxon Rank-Sum Test: 
# Test indicate distributions are not significantly different
#-------COMPARING BIOLEVELS-GCF WITH BIOLEVELS-RANK DIVERSITY FUNCTION --------------

#Plots biome-taxa against biome-gcf-shannon at a certain biome level
plot_biolevel_gcf_rank <- function(level_gcf_diversity, level_genus_diversity, 
                                   level_family_diversity, level_species_diversity, level) {
  
  # Convert level to symbol for tidy evaluation
  level_sym <- ensym(level)
  
  # Merge all dataframes on the selected level column
  merged_data <- level_gcf_diversity %>%
    inner_join(level_genus_diversity, by = rlang::as_name(level_sym), suffix = c("", "_genus")) %>%
    inner_join(level_family_diversity, by = rlang::as_name(level_sym), suffix = c("", "_family")) %>%
    inner_join(level_species_diversity, by = rlang::as_name(level_sym), suffix = c("", "_species")) %>%
    arrange(richness_genus)  # Sorting by increasing richness of genus
  
  # Keep only richness columns
  richness_data <- merged_data %>%
    select(richness, richness_genus, richness_family, richness_species)  
  
  # Compute Pearson correlation
  correlation_matrix <- cor(richness_data, method = "spearman")
  
  print(correlation_matrix[2:4,1])
  # Display as a table
  #kable(correlation_matrix, digits = 3, format = "html", caption = "Pearson Correlation Matrix") %>%
  #kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
  
  # Create the plot
  ggplot(merged_data, aes(x = reorder(!!level_sym, richness_genus))) +
    # Points
    geom_point(aes(y = richness, color = "GCF"), size = 3) + 
    geom_point(aes(y = richness_genus, color = "Genus"), size = 3) + 
    geom_point(aes(y = richness_family, color = "Family"), size = 3) +
    geom_point(aes(y = richness_species, color = "Species"), size = 3) +
    
    # Lines
    geom_line(aes(y = richness, group = 1, color = "GCF"), linetype = "dashed") +
    geom_line(aes(y = richness_genus, group = 1, color = "Genus"), linetype = "solid") +
    geom_line(aes(y = richness_family, group = 1, color = "Family"), linetype = "dotdash") +
    geom_line(aes(y = richness_species, group = 1, color = "Species"), linetype = "twodash") +
    
    # Labels and theme
    labs(x = as_label(level_sym), y = "Richness", color = "Diversity Type") +
    scale_color_manual(values = c("GCF" = "tomato", 
                                  "Genus" = "turquoise3", 
                                  "Family" = "turquoise4", 
                                  "Species" = "turquoise2")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -60, hjust = 0))
  
    
}

#Plots mean-sample-taxa-shannon against biome-gcf-shannon at a certain biome level
plot_sample_gcf_diversity <- function(sample_genus_diversity, sample_species_diversity,
                                      sample_family_diversity, sample_gcf_diversity, level) {
  
  # Combine all datasets into one unified dataframe
  combined_data <- bind_rows(
    sample_genus_diversity %>% mutate(type = "Genus"),
    sample_species_diversity %>% mutate(type = "Species"),
    sample_family_diversity %>% mutate(type = "Family"),
    sample_gcf_diversity %>% mutate(type = "GCF") 
  )
  
  # Compute mean and standard error for each type and biome level
  summary_data <- combined_data %>%
    group_by(type, !!ensym(level)) %>%
    summarise(
      mean_shannon = mean(shannon, na.rm = TRUE),
      se_shannon = sd(shannon, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup()
  
  # Compute the order based on increasing Shannon value for GCF
  ordering <- summary_data %>%
    filter(type == "GCF") %>%
    arrange(mean_shannon) %>%
    pull(!!ensym(level))
  
  # Reorder the x-axis based on Shannon order in GCF
  summary_data <- summary_data %>%
    mutate(!!ensym(level) := factor(!!ensym(level), levels = ordering))
  
  #Remove NA
  summary_data <- summary_data %>% 
    filter(!is.na(!!ensym(level)))
  
  # Create the plot
 plot <- ggplot(summary_data, aes(x = !!ensym(level), y = mean_shannon, color = type)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon),
                  width = 0.2) +
    geom_line(aes(group = type), size = 1, linetype = "solid") +
    
    # Labels and theme
    labs(x = as_label(ensym(level)), y = "Shannon Index", color = "Diversity Type") +
    
    # Custom color palette
    scale_color_manual(values = c("GCF" = "tomato", 
                                  "Genus" = "turquoise3", 
                                  "Family" = "turquoise4", 
                                  "Species" = "turquoise2")) +
    
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -60, hjust = 0))
  
  print(plot)
  
  return(summary_data)
}


#Plot that includes counts V.2
plot_sample_gcf_indexes <- function(genus_diversity, species_diversity,
                                    family_diversity, gcf_diversity, 
                                    gcf_counts, level, order = "d") {
  
  
  # Step 1: Merge datasets progressively using assembly key
  combined_data <- genus_diversity %>%
    inner_join(species_diversity[, c("assembly", "shannon")], by = "assembly", suffix = c("_genus", "_species")) %>%
    inner_join(family_diversity[, c("assembly", "shannon")], by = "assembly", suffix = c("", "_family")) %>%
    inner_join(gcf_diversity[, c("assembly", "shannon")], by = "assembly", suffix = c("", "_gcf")) %>%
    inner_join(gcf_counts[, c("assembly", "counts")], by = "assembly")
  
  # Step 2: Group and calculate mean and standard error
  summary_data <- combined_data %>%
    group_by(!!ensym(level)) %>%
    summarise(
      mean_family_shannon = mean(shannon, na.rm = TRUE),
      se_family_shannon = sd(shannon, na.rm = TRUE) / sqrt(n()),
      
      mean_genus_shannon = mean(shannon_genus, na.rm = TRUE),
      se_genus_shannon = sd(shannon_genus, na.rm = TRUE) / sqrt(n()),
      
      mean_species_shannon = mean(shannon_species, na.rm = TRUE),
      se_species_shannon = sd(shannon_species, na.rm = TRUE) / sqrt(n()),
      
      mean_gcf_shannon = mean(shannon_gcf, na.rm = TRUE),
      se_gcf_shannon = sd(shannon_gcf, na.rm = TRUE) / sqrt(n()),
      
      mean_counts = mean(counts, na.rm = TRUE),
      se_counts = sd(counts, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup()
  
  # Step 3: Handle ordering based on the `order` argument
  if (order == "d") {
    # Order by increasing mean_counts
    summary_data <- summary_data %>%
      arrange(mean_counts) %>%
      mutate(!!ensym(level) := factor(!!ensym(level), levels = !!ensym(level)))
  } else if (order == "a") {
    # Order alphabetically by the `level` column
    summary_data <- summary_data %>%
      arrange(!!ensym(level)) %>%
      mutate(!!ensym(level) := factor(!!ensym(level), levels = unique(!!ensym(level))))
  }
  
  # Step 4: Remove NA values
  summary_data <- summary_data %>% filter(!is.na(!!ensym(level)))
  
  
  # Step 5: Create plot
  p <- ggplot(summary_data, aes(x = !!ensym(level))) +
    # Counts on secondary axis (rescaled for better visibility)
    geom_point(aes(y = mean_counts / 20, color = "GCF Richness"), size = 3, shape = 17) +
    geom_line(aes(y = mean_counts / 20, group = 1, color = "GCF Richness"), size = 1, linetype = "dotted") +
    geom_errorbar(aes(ymin = (mean_counts - se_counts) / 20, ymax = (mean_counts + se_counts) / 20, color = "GCF Richness"),
                  width = 0.2) +
    
    # GCF shannon
    geom_point(aes(y = mean_gcf_shannon, color = "GCF Shannon"), size = 3) +
    geom_line(aes(y = mean_gcf_shannon, group = 1, color = "GCF Shannon"), size = 1) +
    geom_errorbar(aes(ymin = mean_gcf_shannon - se_gcf_shannon, ymax = mean_gcf_shannon + se_gcf_shannon, color = "GCF Shannon"),
                  width = 0.2) +
    
    # Family shannon
    geom_point(aes(y = mean_family_shannon, color = "Family Shannon"), size = 3) +
    geom_line(aes(y = mean_family_shannon, group = 1, color = "Family Shannon"), size = 1) +
    geom_errorbar(aes(ymin = mean_family_shannon - se_family_shannon, ymax = mean_family_shannon + se_family_shannon, color = "Family Shannon"),
                  width = 0.2) +
    
    # Genus shannon
    geom_point(aes(y = mean_genus_shannon, color = "Genus Shannon"), size = 3) +
    geom_line(aes(y = mean_genus_shannon, group = 1, color = "Genus Shannon"), size = 1) +
    geom_errorbar(aes(ymin = mean_genus_shannon - se_genus_shannon, ymax = mean_genus_shannon + se_genus_shannon, color = "Genus Shannon"),
                  width = 0.2) +
    
    # Species shannon
    geom_point(aes(y = mean_species_shannon, color = "Species Shannon"), size = 3) +
    geom_line(aes(y = mean_species_shannon, group = 1, color = "Species Shannon"), size = 1) +
    geom_errorbar(aes(ymin = mean_species_shannon - se_species_shannon, ymax = mean_species_shannon + se_species_shannon, color = "Species Shannon"),
                  width = 0.2) +
    
    # Scale for second y-axis (counts)
    scale_y_continuous(
      name = "Shannon",
      sec.axis = sec_axis(~ . * 20, name = "GCF Richness") # Reverse the rescaling
    ) +
    
    # Labels and theme
    #labs( x = as.character(ensym(level)), title = paste("Mean Shannon and Richness by", as.character(ensym(level)))) +
    labs( x = as.character(ensym(level))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -60, hjust = 0)) +
    
    # Add legend
    scale_color_manual(
      name = "Diversity",
      values = c(
        "GCF Richness" = "gray40",
        "GCF Shannon" = "tomato",
        "Family Shannon" = "turquoise4",
        "Genus Shannon" = "turquoise3",
        "Species Shannon" = "turquoise2"
      )
    )
  
  print(p)
  
  return(summary_data)
  
}


# Function to plot richness and counts by biome level
plot_sample_gcf_richness <- function(sample_taxa_richness, sample_gcf_counts, biome_level) {
  
  # Convert richness columns to numeric (if necessary)
  sample_taxa_richness <- sample_taxa_richness %>%
    mutate(
      family_richness = as.numeric(as.character(family_richness)),
      genus_richness = as.numeric(as.character(genus_richness)),
      species_richness = as.numeric(as.character(species_richness))
    )
  
  # Merge datasets by shared columns
  combined_data <- inner_join(sample_taxa_richness, sample_gcf_counts[, c("assembly", "counts")], by = c("assembly"))
  
  # Group and calculate mean and standard error
  summary_data <- combined_data %>%
    group_by(!!ensym(biome_level)) %>%
    summarise(
      mean_family_richness = mean(family_richness, na.rm = TRUE),
      se_family_richness = sd(family_richness, na.rm = TRUE) / sqrt(n()),
      
      mean_genus_richness = mean(genus_richness, na.rm = TRUE),
      se_genus_richness = sd(genus_richness, na.rm = TRUE) / sqrt(n()),
      
      mean_species_richness = mean(species_richness, na.rm = TRUE),
      se_species_richness = sd(species_richness, na.rm = TRUE) / sqrt(n()),
      
      mean_counts = mean(counts, na.rm = TRUE),
      se_counts = sd(counts, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup() %>%
    arrange(mean_counts) # Order by increasing counts
  
  #Remove NA
  summary_data <- summary_data %>% filter(!is.na(!!ensym(biome_level)))
  
  # Create plot
  p <- ggplot(summary_data, aes(x = reorder(!!ensym(biome_level), mean_counts))) +
    # Counts on secondary axis (rescaled for better visibility)
    geom_point(aes(y = mean_counts * 50, color = "GCF Richness"), size = 3) +
    geom_line(aes(y = mean_counts * 50, group = 1, color = "GCF Richness"), size = 1) +
    geom_errorbar(aes(ymin = (mean_counts - se_counts) * 50, ymax = (mean_counts + se_counts) * 50, color = "GCF Richness"),
                  width = 0.2) +
    
    # Family richness
    geom_point(aes(y = mean_family_richness, color = "Family Richness"), size = 3) +
    geom_line(aes(y = mean_family_richness, group = 1, color = "Family Richness"), size = 1) +
    geom_errorbar(aes(ymin = mean_family_richness - se_family_richness, ymax = mean_family_richness + se_family_richness, color = "Family Richness"),
                  width = 0.2) +
    
    # Genus richness
    geom_point(aes(y = mean_genus_richness, color = "Genus Richness"), size = 3) +
    geom_line(aes(y = mean_genus_richness, group = 1, color = "Genus Richness"), size = 1) +
    geom_errorbar(aes(ymin = mean_genus_richness - se_genus_richness, ymax = mean_genus_richness + se_genus_richness, color = "Genus Richness"),
                  width = 0.2) +
    
    # Species richness
    geom_point(aes(y = mean_species_richness, color = "Species Richness"), size = 3) +
    geom_line(aes(y = mean_species_richness, group = 1, color = "Species Richness"), size = 1) +
    geom_errorbar(aes(ymin = mean_species_richness - se_species_richness, ymax = mean_species_richness + se_species_richness, color = "Species Richness"),
                  width = 0.2) +
    
    # Scale for second y-axis (counts)
    scale_y_continuous(
      name = "Mean Taxa Richness",
      sec.axis = sec_axis(~ . / 50, name = "Mean GCF Richness") # Reverse the rescaling
    ) +
    
    # Labels and theme
    #labs(x = as.character(ensym(biome_level)), title = paste("Mean Richness by", as.character(ensym(biome_level))) ) +
    labs(x = as.character(ensym(biome_level)) ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -60, hjust = 0)) +
    
    # Add legend
    scale_color_manual(
      name = "Diversity",
      values = c(
        "GCF Richness" = "red3",
        "Family Richness" = "darkolivegreen4",
        "Genus Richness" = "darkolivegreen3",
        "Species Richness" = "darkolivegreen2"
      )
    )
  
  print(p)
  
  return(summary_data)
}


#

#---------COMPARING TAXA-DIVERSITY AGAINST GCF-DIVERSITY (REGIONS)--------

#Comparing using only region biome-taxa diversity
plot_biolevel_gcf_rank(biome_gcf_diversity, biome_genus_diversity, biome_fam_diversity, biome_spec_diversity, biome)

plot_biolevel_gcf_rank(subbiome_gcf_diversity[!(is.na(subbiome_gcf_diversity$subbiome)), ], subbiome_genus_diversity[!(is.na(subbiome_genus_diversity$subbiome)), ]
                       , subbiome_fam_diversity[!(is.na(subbiome_fam_diversity$subbiome)), ], subbiome_spec_diversity[!(is.na(subbiome_spec_diversity$subbiome)), ], subbiome)

plot_biolevel_gcf_rank(infrabiome_gcf_diversity[!(is.na(infrabiome_gcf_diversity$infrabiome)), ], infrabiome_genus_diversity[!(is.na(infrabiome_genus_diversity$infrabiome)), ]
                       , infrabiome_fam_diversity[!(is.na(infrabiome_fam_diversity$infrabiome)), ], infrabiome_spec_diversity[!(is.na(infrabiome_spec_diversity$infrabiome)), ], infrabiome)

plot_biolevel_gcf_rank(microbiome_gcf_diversity[!(is.na(microbiome_gcf_diversity$microbiome)), ], microbiome_genus_diversity[!(is.na(microbiome_genus_diversity$microbiome)), ]
                       , microbiome_fam_diversity[!(is.na(microbiome_fam_diversity$microbiome)), ], microbiome_spec_diversity[!(is.na(microbiome_spec_diversity$microbiome)), ], microbiome)

#---------COMPARING TAXA-SHANNON AGAINST GCF-SHANNON (SAMPLE MEAN)--------

#Comparing using the full set sample-taxa shannon 
richness_biome_summary <- plot_sample_gcf_diversity(sample_genus_diversity[!(is.na(sample_genus_diversity$biome)), ], sample_species_diversity[!(is.na(sample_species_diversity$biome)), ], 
                          sample_family_diversity[!(is.na(sample_family_diversity$biome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$biome)), ], biome)

summary_subbiome <- plot_sample_gcf_diversity(sample_genus_diversity[!(is.na(sample_genus_diversity$subbiome)), ], sample_species_diversity[!(is.na(sample_species_diversity$subbiome)), ],
                          sample_family_diversity[!(is.na(sample_family_diversity$subbiome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$subbiome)), ], subbiome) 

summary_infrabiome <- plot_sample_gcf_diversity(sample_genus_diversity[!(is.na(sample_family_diversity$infrabiome)), ], sample_species_diversity[!(is.na(sample_species_diversity$infrabiome)), ], 
                          sample_family_diversity[!(is.na(sample_family_diversity$infrabiome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$infrabiome)), ], infrabiome) 

summary_microbiome <- plot_sample_gcf_diversity(sample_genus_diversity[!(is.na(sample_family_diversity$microbiome)), ], sample_species_diversity[!(is.na(sample_species_diversity$microbiome)), ], 
                          sample_family_diversity[!(is.na(sample_family_diversity$microbiome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$microbiome)), ], microbiome) 
#

#---------COMPARING TAXA-SHANNON, GCF-SHANNON, GCF-RICHNESS (SAMPLE MEAN)---------


shannon_biome_summary <- plot_sample_gcf_indexes(sample_genus_diversity[!(is.na(sample_genus_diversity$biome)), ], sample_species_diversity[!(is.na(sample_species_diversity$biome)), ], 
                                           sample_family_diversity[!(is.na(sample_family_diversity$biome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$biome)), ],
                                           sample_gcf_counts[!(is.na(sample_gcf_counts$biome)), ], biome, order = "d")

shannon_subbiome_summary <- plot_sample_gcf_indexes(sample_genus_diversity[!(is.na(sample_genus_diversity$subbiome)), ], sample_species_diversity[!(is.na(sample_species_diversity$subbiome)), ], 
                        sample_family_diversity[!(is.na(sample_family_diversity$subbiome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$subbiome)), ],
                        sample_gcf_counts[!(is.na(sample_gcf_counts$subbiome)), ], subbiome, order = "d")

shannon_infrabiome_summary <- plot_sample_gcf_indexes(sample_genus_diversity[!(is.na(sample_genus_diversity$infrabiome)), ], sample_species_diversity[!(is.na(sample_species_diversity$infrabiome)), ], 
                        sample_family_diversity[!(is.na(sample_family_diversity$infrabiome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$infrabiome)), ],
                        sample_gcf_counts[!(is.na(sample_gcf_counts$infrabiome)), ], infrabiome, order = "d")

shannon_microbiome_summary <- plot_sample_gcf_indexes(sample_genus_diversity[!(is.na(sample_genus_diversity$microbiome)), ], sample_species_diversity[!(is.na(sample_species_diversity$microbiome)), ], 
                        sample_family_diversity[!(is.na(sample_family_diversity$microbiome)), ], sample_gcf_diversity[!(is.na(sample_gcf_diversity$microbiome)), ],
                        sample_gcf_counts[!(is.na(sample_gcf_counts$microbiome)), ], microbiome)
#

#---------COMPARING TAXA-RICHNESS, GCF-RICHNESS (SAMPLE MEAN)------------

richness_biome_summary <- plot_sample_gcf_richness(sample_taxa_richness, sample_gcf_counts, biome)
richness_subbiome_summary <- plot_sample_gcf_richness(sample_taxa_richness, sample_gcf_counts, subbiome)
richness_infrabiome_summary <- plot_sample_gcf_richness(sample_taxa_richness, sample_gcf_counts, infrabiome)
richness_microbiome_summary <- plot_sample_gcf_richness(sample_taxa_richness, sample_gcf_counts, microbiome)

#

#-------CORRELATION FUNCTIONS---------
test_genus_correlation <- function(summary, level) {
  # Filter only GCF and Genus
  test_filtered <- summary %>%
    filter(type %in% c("GCF", "Genus")) %>%
    select(type, !!ensym(level), mean_shannon)%>%
    filter(complete.cases(type, !!ensym(level)))  # Remove missing values
  
  # Pivot data to wide format
  test_wide <- test_filtered %>%
    pivot_wider(names_from = type, values_from = mean_shannon) 
  
  # Perform correlation test
    print(cor.test(test_wide$GCF, test_wide$Genus, method = "pearson"))
    print(cor.test(test_wide$GCF, test_wide$Genus, method = "spearman"))

}

test_family_correlation <- function(summary, level) {
  # Filter only GCF and Genus
  test_filtered <- summary %>%
    filter(type %in% c("GCF", "Family")) %>%
    select(type, !!ensym(level), mean_shannon)%>%
    filter(complete.cases(type, !!ensym(level)))  # Remove missing values
  
  # Pivot data to wide format
  test_wide <- test_filtered %>%
    pivot_wider(names_from = type, values_from = mean_shannon) %>%
    filter(complete.cases(GCF, Family))  # Remove missing values
  
  # Perform correlation test
  print(cor_test <- cor.test(test_wide$GCF, test_wide$Family, method = "pearson"))
  print(cor_test <- cor.test(test_wide$GCF, test_wide$Family, method = "spearman"))
}

test_species_correlation <- function(summary, level) {
  # Filter only GCF and Genus
  test_filtered <- summary %>%
    filter(type %in% c("GCF", "Species")) %>%
    select(type, !!ensym(level), mean_shannon)%>%
    filter(complete.cases(type, !!ensym(level)))  # Remove missing values
  
  # Pivot data to wide format
  test_wide <- test_filtered %>%
    pivot_wider(names_from = type, values_from = mean_shannon) %>%
    filter(complete.cases(GCF, Species))  # Remove missing values
  
  # Perform correlation test
  print(cor_test <- cor.test(test_wide$GCF, test_wide$Species, method = "pearson"))
  print(cor_test <- cor.test(test_wide$GCF, test_wide$Species, method = "spearman"))
}

test_correlation <- function(summary) {

  print("Species")
  # Perform correlation test
  print(cor_test <- cor.test(summary$mean_counts, summary$mean_species_richness, method = "pearson"))
  print(cor_test <- cor.test(summary$mean_counts, summary$mean_species_richness, method = "spearman"))
  
  print("Genus")
  # Perform correlation test
  print(cor_test <- cor.test(summary$mean_counts, summary$mean_genus_richness, method = "pearson"))
  print(cor_test <- cor.test(summary$mean_counts, summary$mean_genus_richness, method = "spearman"))
  
  print("Family")
  # Perform correlation test
  print(cor_test <- cor.test(summary$mean_counts, summary$mean_family_richness, method = "pearson"))
  print(cor_test <- cor.test(summary$mean_counts, summary$mean_family_richness, method = "spearman"))
}

#---------CORRELATION MEAN SHANNON-------
test_species_correlation(summary_biome,  biome)
test_genus_correlation(summary_biome,  biome)
test_family_correlation(summary_biome,  biome)

test <-summary_biome %>% filter(!.[[2]] %in% c("Amphibia", "Mollusca", "Cnidaria"))
test_genus_correlation(test,  biome)

test_species_correlation(summary_subbiome,  subbiome)
test_genus_correlation(summary_subbiome,  subbiome)
test_family_correlation(summary_subbiome,  subbiome)

test_species_correlation(summary_infrabiome,  infrabiome) 
test_genus_correlation(summary_infrabiome,  infrabiome) 
test_family_correlation(summary_infrabiome,  infrabiome) 

test_species_correlation(summary_microbiome,  microbiome)
test_genus_correlation(summary_microbiome,  microbiome) 
test_family_correlation(summary_microbiome,  microbiome)

#---------CORRELATION MEAN RICHNESS-------
test_correlation(richness_biome_summary)
test_correlation(richness_subbiome_summary)
test_correlation(richness_infrabiome_summary)
test_correlation(richness_microbiome_summary)


#-------SUPERKINGDOM SPLIT FUNCTION--------
plot_superbiome_shannon <- function(data, subbiome_col = "subbiome", 
                                    genus_col = "mean_genus_shannon", 
                                    gcf_col = "mean_gcf_shannon") {
  
  # Ensure column names are properly referenced
  subbiome_sym <- rlang::sym(subbiome_col)
  genus_sym <- rlang::sym(genus_col)
  gcf_sym <- rlang::sym(gcf_col)
  
  # Order the data by increasing GCF Shannon index
  data_ordered <- data %>%
    arrange(!!gcf_sym) %>%
    mutate(!!subbiome_sym := factor(!!subbiome_sym, levels = !!subbiome_sym))
  
  # Pivot to long format
  data_long <- data_ordered %>%
    pivot_longer(cols = c(!!genus_sym, !!gcf_sym),
                 names_to = "type", values_to = "shannon_value")
  
  # Plot
  ggplot(data_long, aes(y = shannon_value, x = !!subbiome_sym, color = type, group = type)) +
    geom_line() +
    geom_point() +
    labs(y = "Shannon Index", x = "Subbiome", color = "Shannon Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -60, hjust = 0))
}


#---------CORRELATION MEAN RICHNNES AND SHANNON AT SUPERBIOME--------

env <- c( "Aquatic:Aquaculture", "Aquatic:Estuary" , "Aquatic:Freshwater","Aquatic:Lentic",
          "Aquatic:Marine" ,"Aquatic:Non-marine Saline and Alkaline", "Aquatic:Sediment",
          "Aquatic:Thermal springs" ,"Terrestrial:Soil")

ha <- c("Algae:Red algae", "Animal:Digestive system", "Arthropoda:Digestive system", "Birds:Digestive system", 
        "Fish:Digestive system", "Human:Digestive system", "Human:Skin", "Insecta:Digestive system", 
        "Mammals:Digestive system" , "Mammals:Gastrointestinal tract", "Mammals:Respiratory system", 
        "Microbial:Bacteria", "Plants:Phylloplane", "Plants:Rhizosphere", "Plants:Root")

#Subbiome subset for environmental
richness_subbiome_env_summary <- richness_subbiome_summary[richness_subbiome_summary$subbiome %in% env , ]
shannon_subbiome_env_summary <- shannon_subbiome_summary[shannon_subbiome_summary$subbiome %in% env , ]

print(cor.test(richness_subbiome_env_summary$mean_genus_richness, richness_subbiome_env_summary$mean_counts, method = "pearson"))
#t = 2.418, df = 7, p-value = 0.04623, cor = 0.6746138
print(cor.test(shannon_subbiome_env_summary$mean_genus_shannon, shannon_subbiome_env_summary$mean_gcf_shannon, method = "pearson"))
#t = 2.2387, df = 7, p-value = 0.06019, cor = 0.645939 

plot_superbiome_shannon(richness_subbiome_env_summary, genus_col = "mean_genus_richness", gcf_col = "mean_counts")
plot_superbiome_shannon(shannon_subbiome_env_summary)


#Subbiome subset for host associated
richness_subbiome_ha_summary <- richness_subbiome_summary[richness_subbiome_summary$subbiome %in% ha , ]
shannon_subbiome_ha_summary <- shannon_subbiome_summary[shannon_subbiome_summary$subbiome %in% ha , ]

print(cor.test(richness_subbiome_ha_summary$mean_genus_richness, richness_subbiome_ha_summary$mean_counts, method = "pearson"))
#t = 3.0596, df = 13, p-value = 0.00913,  cor = 0.6470202
print(cor.test(shannon_subbiome_ha_summary$mean_genus_shannon, shannon_subbiome_ha_summary$mean_gcf_shannon, method = "pearson"))
#t = 2.215, df = 13, p-value = 0.04524, cor = 0.5234406 

plot_superbiome_shannon(richness_subbiome_ha_summary, genus_col = "mean_genus_richness", gcf_col = "mean_counts")
plot_superbiome_shannon(shannon_subbiome_ha_summary)

#Comparing distributions
boxplot(shannon_subbiome_env_summary$mean_gcf_shannon, shannon_subbiome_ha_summary$mean_gcf_shannon)
ks.test(shannon_subbiome_env_summary$mean_gcf_shannon, shannon_subbiome_ha_summary$mean_gcf_shannon)
#D = 0.22222, p-value = 0.8987
wilcox.test(shannon_subbiome_env_summary$mean_gcf_shannon, shannon_subbiome_ha_summary$mean_gcf_shannon)
#W = 62, p-value = 0.7702
