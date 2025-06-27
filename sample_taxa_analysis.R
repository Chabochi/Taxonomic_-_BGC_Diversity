#------PREPROCESSING FUNCTIONS-------

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

#------DATA PROCESSING-----------

#Read assembly to biome
assembly2biome <- read.delim(file = "atlas_tables/assembly2longestbiome.tsv", sep = "\t")

#Read genus diversity of each sample
sample_genus_diversity <- read.delim(file = "atlas_tables/alpha_diversity_results.csv", sep = ",", header=FALSE)
colnames(sample_genus_diversity) <- c("assembly", "shannon")
# Perform the join and add the new column
sample_genus_diversity <- sample_genus_diversity %>%
  left_join(assembly2biome %>% select(assembly, longest_biome), by = "assembly")
#Adjusting biome
sample_genus_diversity <- adjust_biome(sample_genus_diversity)
#Filter out Mixed, NA, and Engineered biomes
sample_genus_diversity <- sample_genus_diversity[sample_genus_diversity$superbiome %in% c("Host-associated", "Environmental"), ]
#Before:34,836 , After: 31,575, Lost: 9%


#Read species diversity of each sample
sample_species_diversity <- read.delim(file = "atlas_tables/species_diversity_results.csv", sep = ",", header=FALSE)
colnames(sample_species_diversity) <- c("assembly", "shannon")
# Perform the join and add the new column
sample_species_diversity <- sample_species_diversity %>%
  left_join(assembly2biome %>% select(assembly, longest_biome), by = "assembly")
#Adjusting biome
sample_species_diversity <- adjust_biome(sample_species_diversity)
#Filter out Mixed, NA, and Engineered biomes
sample_species_diversity <- sample_species_diversity[sample_species_diversity$superbiome %in% c("Host-associated", "Environmental"), ]
#Before:34,836 , After: 31,575, Lost: 9%


#Read family diversity of each sample
sample_family_diversity <- read.delim(file = "atlas_tables/family_diversity_results.csv", sep = ",", header=FALSE)
colnames(sample_family_diversity) <- c("assembly", "shannon")
# Perform the join and add the new column
sample_family_diversity <- sample_family_diversity %>%
  left_join(assembly2biome %>% select(assembly, longest_biome), by = "assembly")
#Adjusting biome
sample_family_diversity <- adjust_biome(sample_family_diversity)
#Filter out Mixed, NA, and Engineered biomes
sample_family_diversity <- sample_family_diversity[sample_family_diversity$superbiome %in% c("Host-associated", "Environmental"), ]
#Before:34,836 , After: 31,577, Lost: 9%

#Read  gcf shannon index of each sample
sample_gcf_diversity <- read.delim(file = "atlas_tables/sample_gcf_shannon.tsv", sep = "\t", header=FALSE)
colnames(sample_gcf_diversity) <- c("assembly", "shannon")
# Perform the join and add the new column
sample_gcf_diversity <- sample_gcf_diversity %>%
  left_join(assembly2biome %>% select(assembly, longest_biome), by = "assembly")
#Adjusting biome
sample_gcf_diversity <- adjust_biome(sample_gcf_diversity)
sample_gcf_diversity <- sample_gcf_diversity[sample_gcf_diversity$superbiome %in% c("Host-associated", "Environmental"), ]
#Before:29,632 , After:, Lost:26,641 

#Read adjusted richness for gcf diversity of each sample
sample_gcf_counts <- read.delim(file = "atlas_tables/sample_gcf_counts.tsv", sep = "\t", header=FALSE)
colnames(sample_gcf_counts) <- c("assembly", "counts")
# Perform the join and add the new column
sample_gcf_counts <- sample_gcf_counts %>%
  left_join(assembly2biome %>% select(assembly, longest_biome), by = "assembly")
#Adjusting biome
sample_gcf_counts <- adjust_biome(sample_gcf_counts)
sample_gcf_counts <- sample_gcf_counts[sample_gcf_counts$superbiome %in% c("Host-associated", "Environmental"), ]
#Before:29,627 , After:, Lost:26,636


#Read richness for taxa 
sample_taxa_richness <- read.delim(file = "atlas_tables/sample_taxa_richness.tsv", sep = "\t", header=TRUE)
# Perform the join and add the new column
sample_taxa_richness <- sample_taxa_richness %>%
  left_join(assembly2biome %>% select(assembly, longest_biome), by = "assembly")
#Adjusting biome
sample_taxa_richness <- adjust_biome(sample_taxa_richness)
sample_taxa_richness <- sample_taxa_richness[sample_taxa_richness$superbiome %in% c("Host-associated", "Environmental"), ]
#Before:34,841 , After:, Lost:31,581

#------SUMMARIZING FUNCTIONS-----

summarize_sample <- function(data, level, rank = "g"){
  level_sym <- ensym(level)
  
  if (rank == "s") {
    data_diversity_mean <- data %>%
      group_by(!!level_sym) %>%
      summarise(species_mean = mean(shannon),
                species_std = sd(shannon))
  } else if (rank == "g") {
    data_diversity_mean <- data %>%
      group_by(!!level_sym) %>%
      summarise(genus_mean = mean(shannon),
                genus_std = sd(shannon))
  } else if (rank == "f") {
  data_diversity_mean <- data %>%
    group_by(!!level_sym) %>%
    summarise(family_mean = mean(shannon),
              family_std = sd(shannon))
  } else if (rank == "") {
    data_diversity_mean <- data %>%
      group_by(!!level_sym) %>%
      summarise(mean = mean(shannon),
                std = sd(shannon))}
  
  return(data_diversity_mean)
}

#IN PROCESS
merge_mean_sample <- function(genus, family, species, level){
  level_sym <- ensym(level)
  
  sample_genus_diversity_mean <- summarize_sample(genus, level)
  sample_family_diversity_mean <- summarize_sample(family, level, "f") 
  sample_species_diversity_mean <- summarize_sample(species, level, "s")
  
  #merge all data frames in list
  merge_data <- list(sample_genus_diversity_mean, sample_family_diversity_mean, sample_species_diversity_mean) %>% 
    reduce(full_join, by=as.character(level))
  
  return(merge_data)
}


#---------SUMMARIZING BY BIOME------
#Summarising sample taxa diversity
sample_genus_diversity_mean <- summarize_sample(sample_genus_diversity, biome)
sample_family_diversity_mean <- summarize_sample(sample_family_diversity, biome, "f") 
sample_species_diversity_mean <- summarize_sample(sample_species_diversity, biome, "s")

#merge all data frames in list
biome_sample_taxa_diversity <- list(sample_genus_diversity_mean, sample_family_diversity_mean, sample_species_diversity_mean) %>% 
  reduce(full_join, by='biome')

#Summarising sample gcf diversity
sample_biome_gcf_diversity_mean <- summarize_sample(sample_gcf_diversity, biome, "")
sample_subbiome_gcf_diversity_mean <- summarize_sample(sample_gcf_diversity, subbiome, "")
sample_infrabiome_gcf_diversity_mean <- summarize_sample(sample_gcf_diversity, infrabiome, "")
sample_microbiome_gcf_diversity_mean <- summarize_sample(sample_gcf_diversity, microbiome, "")
#

#---------SUMMARIZING BY SUBBIOME------
sample_genus_diversity_mean <- summarize_sample(sample_genus_diversity, subbiome)
sample_family_diversity_mean <- summarize_sample(sample_family_diversity, subbiome, "f") 
sample_species_diversity_mean <- summarize_sample(sample_species_diversity, subbiome, "s")

#merge all data frames in list
subbiome_sample_taxa_diversity <- list(sample_genus_diversity_mean, sample_family_diversity_mean, sample_species_diversity_mean) %>% 
  reduce(full_join, by='subbiome')

#
#---------SUMMARIZING BY INFRABIOME------
sample_genus_diversity_mean <- summarize_sample(sample_genus_diversity, infrabiome)
sample_family_diversity_mean <- summarize_sample(sample_family_diversity, infrabiome, "f") 
sample_species_diversity_mean <- summarize_sample(sample_species_diversity, infrabiome, "s")

#merge all data frames in list
infrabiome_sample_taxa_diversity <- list(sample_genus_diversity_mean, sample_family_diversity_mean, sample_species_diversity_mean) %>% 
  reduce(full_join, by='infrabiome')
#

#---------SUMMARIZING BY MICROBIOME------
sample_genus_diversity_mean <- summarize_sample(sample_genus_diversity, microbiome)
sample_family_diversity_mean <- summarize_sample(sample_family_diversity, microbiome, "f") 
sample_species_diversity_mean <- summarize_sample(sample_species_diversity, microbiome, "s")

#merge all data frames in list
microbiome_sample_taxa_diversity <- list(sample_genus_diversity_mean, sample_family_diversity_mean, sample_species_diversity_mean) %>% 
  reduce(full_join, by='microbiome')
#

#---------COMPARING BIOMELEVEL MEAN-SAMPLE-GENUS SHANNON-------

# Create a combined dataframe
biome_level_sample_genus_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_sample_taxa_diversity$genus_mean),
  data.frame(Level = "Subbiome", Shannon = subbiome_sample_taxa_diversity$genus_mean),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_sample_taxa_diversity$genus_mean),
  data.frame(Level = "Microbiome", Shannon = microbiome_sample_taxa_diversity$genus_mean)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order


# Create the boxplot
ggplot(biome_level_sample_genus_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Genus Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_sample_taxa_diversity$genus_mean,microbiome_sample_taxa_diversity$genus_mean)
# Kolmogorov-Smirnov Test: D = 0.27843, p-value = 0.3217
wilcox.test(biome_sample_taxa_diversity$genus_mean,microbiome_sample_taxa_diversity$genus_mean)
# Wilcoxon Rank-Sum Test: W = 218, p-value = 0.4327
# Test indicate distributions are not significantly different

#---------COMPARING BIOMELEVEL MEAN-SAMPLE-SPECIES SHANNON-------
# Create a combined dataframe
biome_level_sample_species_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_sample_taxa_diversity$species_mean),
  data.frame(Level = "Subbiome", Shannon = subbiome_sample_taxa_diversity$species_mean),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_sample_taxa_diversity$species_mean),
  data.frame(Level = "Microbiome", Shannon = microbiome_sample_taxa_diversity$species_mean)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

# Create the boxplot
ggplot(biome_level_sample_species_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Species Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_sample_taxa_diversity$species_mean,microbiome_sample_taxa_diversity$species_mean)
# Kolmogorov-Smirnov Test: D = 0.16078, p-value = 0.8995
wilcox.test(biome_sample_taxa_diversity$species_mean,microbiome_sample_taxa_diversity$species_mean)
# Wilcoxon Rank-Sum Test: W = 243, p-value = 0.8052
# Test indicate distributions are not significantly different

#---------COMPARING BIOMELEVEL MEAN-SAMPLE-FAMILY SHANNON-------
# Create a combined dataframe
biome_level_sample_family_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_sample_taxa_diversity$family_mean),
  data.frame(Level = "Subbiome", Shannon = subbiome_sample_taxa_diversity$family_mean),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_sample_taxa_diversity$family_mean),
  data.frame(Level = "Microbiome", Shannon = microbiome_sample_taxa_diversity$family_mean)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

# Create the boxplot
ggplot(biome_level_sample_family_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Family Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_sample_taxa_diversity$family_mean,microbiome_sample_taxa_diversity$family_mean)
# Kolmogorov-Smirnov Test: D = 0.30784, p-value = 0.2199
wilcox.test(biome_sample_taxa_diversity$family_mean,microbiome_sample_taxa_diversity$family_mean)
# Wilcoxon Rank-Sum Test: W = 228, p-value = 0.5694
# Test indicate distributions are not significantly different

#--------COMPARING BIOMELEVEL MEAN-SAMPLE-GCF SHANNON--------
biome_level_sample_gcf_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = sample_biome_gcf_diversity_mean$mean),
  data.frame(Level = "Subbiome", Shannon = sample_subbiome_gcf_diversity_mean$mean),
  data.frame(Level = "Infrabiome", Shannon = sample_infrabiome_gcf_diversity_mean$mean),
  data.frame(Level = "Microbiome", Shannon = sample_microbiome_gcf_diversity_mean$mean)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

# Create the boxplot
ggplot(biome_level_sample_gcf_shannon, aes(x = Level, y = Shannon, fill = "Red1")) +
  geom_boxplot() +
  #scale_fill_brewer(palette = "Reds") +  # Try "Set1", "Set2", "Pastel1", "Dark2", et
  #scale_fill_brewer(palette = "Reds") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "GCF Shannon Index") +
  #axis.text.x = element_text(size=14, angle=45) +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(sample_biome_gcf_diversity_mean$mean,sample_subbiome_gcf_diversity_mean$mean)
# Kolmogorov-Smirnov Test: D = 0.1675, p-value = 0.887
wilcox.test(sample_biome_gcf_diversity_mean$mean,sample_subbiome_gcf_diversity_mean$mean)
# Wilcoxon Rank-Sum Test: W = 182.5, p-value = 0.6495
# Test indicate distributions are not significantly different

#Analyzing if distributions are significantly different
ks.test(sample_biome_gcf_diversity_mean$mean,sample_infrabiome_gcf_diversity_mean$mean)
# Kolmogorov-Smirnov Test: D = 0.32048, p-value = 0.141
wilcox.test(sample_biome_gcf_diversity_mean$mean,sample_infrabiome_gcf_diversity_mean$mean)
# Wilcoxon Rank-Sum Test: W = 306.5, p-value = 0.2759
# Test indicate distributions are not significantly different

#Analyzing if distributions are significantly different
ks.test(sample_biome_gcf_diversity_mean$mean,sample_microbiome_gcf_diversity_mean$mean)
# Kolmogorov-Smirnov Test: D = 0.42235, p-value = 0.02554
wilcox.test(sample_biome_gcf_diversity_mean$mean,sample_microbiome_gcf_diversity_mean$mean)
# Wilcoxon Rank-Sum Test: W = 157, p-value = 0.02203
# Test indicate distributions are not significantly different

#Analyzing if distributions are significantly different
ks.test(sample_infrabiome_gcf_diversity_mean$mean,sample_microbiome_gcf_diversity_mean$mean)
# Kolmogorov-Smirnov Test: D = 0.42235, p-value = 0.02554
wilcox.test(sample_infrabiome_gcf_diversity_mean$mean,sample_microbiome_gcf_diversity_mean$mean)
# Wilcoxon Rank-Sum Test: W = 157, p-value = 0.02203
# Test indicate distributions are not significantly different


#------MULTIRANK FOR MEAN-SAMPLE-TAXA-SHANNON FUNCTION--------------

#Comparing sample-taxa-diversity split by biome
plot_biome_level_shannon <- function(sample_genus_diversity, sample_species_diversity, sample_family_diversity) {
  
  # Combine the dataframes with a new 'type' column
  combined_data <- bind_rows(
    sample_genus_diversity %>% mutate(type = "Genus"),
    sample_species_diversity %>% mutate(type = "Species"),
    sample_family_diversity %>% mutate(type = "Family")
  )
  
  # Average Shannon values for each unique biome at each level
  summary_data <- combined_data %>%
    pivot_longer(cols = c(biome, subbiome, infrabiome, microbiome),
                 names_to = "biome_level", values_to = "biome_name") %>%
    drop_na(biome_name) %>%  # Remove NA values
    group_by(type, biome_level, biome_name) %>%
    summarise(mean_shannon = mean(shannon, na.rm = TRUE), .groups = "drop")
  
  # Reorder the biome levels to reflect hierarchy
  summary_data <- summary_data %>%
    mutate(biome_level = factor(biome_level, levels = c("biome", "subbiome", "infrabiome", "microbiome")))
  
  # Create the box plot
  p <- ggplot(summary_data, aes(x = biome_level, y = mean_shannon, fill = type)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), width = 0.6) +
    labs(x = "Biome Level", y = "Mean Shannon Index", fill = "Taxonomic Rank") +
    
    # Custom color palette
    scale_fill_manual(values = c("Family" = "darkseagreen4", "Genus" = "darkseagreen", "Species" = "darkseagreen2")) +
    
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot(p)
  
}

#---------MULTIRANK FOR MEAN-SAMPLE-TAXA-SHANNON--------------

#Comparing MEAN-SAMPLE-TAXA-SHANNON by rank and biome
plot_biome_level_shannon(sample_genus_diversity, sample_species_diversity, sample_family_diversity)
