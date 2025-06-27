#FINAL ANALYSIS
print("welcome")

#-------LIBRARIES-------
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse) #To parse the biome names
library(forcats)  # Needed for fct_reorder()
library(gridExtra) #to plot multiple graphs
library(stringr) #to modify categories
#library(indicspecies) # For species specificity
#library(data.table) #to calculate specificity of bgc

#------PREPROCESSING FUNCTIONS-------

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

#Used as reference for caluclating biome-gcf-diversity
regions_taxonomy_biome <- adjust_biome(regions_taxonomy)
regions_taxonomy_biome  <- regions_taxonomy_biome [regions_taxonomy_biome$superbiome %in% c("Host-associated", "Environmental"), ]

#Read genus taxdump file
gen_taxdump<-read.delim(file = "atlas_tables/genus_taxdump.tsv", sep = "\t", header=FALSE)
colnames(gen_taxdump) <- c("tax_id","parent","rank","name", "genus_id")

#Creating regions table with genus taxonomy
regions_with_genus_named <-  name_regions(regions_taxonomy, gen_taxdump)

#PROCESSING GCF IDs
#Turning gcf ids from integer to character
regions_with_genus_named[,14] <- sapply(regions_with_genus_named[,14], as.character)

#PROCESSING LONGEST_BIOME#
regions_with_genus_named <- adjust_biome(regions_with_genus_named)
#Filter out Mixed, NA, and Engineered biomes
regions_with_genus_named <- regions_with_genus_named[regions_with_genus_named$superbiome %in% c("Host-associated", "Environmental"), ]
#Before: 1,634,728 , After: 1,457,665, Lost: 10.9 %

#PROCESSING CATEGORY AND TYPE LABELS
# Usage for product_categories
regions_with_genus_named <- adjust_labels(regions_with_genus_named, bigslice_gcf_id, product_categories)
# Usage for type
regions_with_genus_named <- adjust_labels(regions_with_genus_named, bigslice_gcf_id, type)


#SUMMARY
cat("Number of unique genus:", n_distinct(regions_with_genus_named$genus_name), "\n\n",
  "Number of unique GCFs:", n_distinct(regions_with_genus_named$bigslice_gcf_id), "\n\n",
  "Number of unique types:", n_distinct(regions_with_genus_named$type), "\n\n",
  "Number of unique biomes:", n_distinct(regions_with_genus_named$longest_biome), "\n")

#Number of unique genus (2,274)
#Number of unique GCFs (16,041)
#Number of unique types (575)
#Number of unique biomes (97)
#Number of unique assemblies (26,308)

#Remove non used tables
remove(regions_taxonomy, gen_taxdump)

#Write table
#write.table(select(regions_with_genus_named, c("region_id","genus_name","bigslice_gcf_id")), file='region_with_genus.csv', quote=FALSE, sep=',', col.names = FALSE, row.names = FALSE)


assembly2biome_split <- adjust_biome(assembly2biome)
#Superbiome : 35400
#Biome : 33199 (93%)
#Subbiome : 32080 (90%)
#Infrabiome : 21401 (60%)
#Microbiome : 16332 (46%)
remove(assembly2biome_split)



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

#--------BIOME LEVEL COUNTS---------

count_biome_levels <- function(df) {
  # Convert to long format while preserving the biome
  df_long <- df %>%
    pivot_longer(cols = c(biome, subbiome, infrabiome, microbiome), 
                 names_to = "level", values_to = "level_name") %>%
    drop_na(level_name)  # Remove NA values
  
  # Keep only the lowest available level for each region
  df_filtered <- df_long %>%
    group_by(region_id) %>%
    filter(row_number() == n()) %>%
    ungroup()
  
  # Add back the original biome as a separate column
  df_filtered <- df_filtered %>%
    left_join(df %>% select(region_id, biome), by = "region_id")
  
  # Count the number of regions per biome
  biome_counts <- df_filtered %>%
    count(biome, name = "count")  # Count occurrences of each biome
  
  # Merge count data back into df_filtered
  df_filtered <- df_filtered %>%
    left_join(biome_counts, by = "biome")
  
  # Ensure correct stacking order (biome at the bottom, microbiome at the top)
  df_filtered <- df_filtered %>%
    mutate(level = factor(level, levels = c("microbiome", "infrabiome", "subbiome", "biome")),
           biome = reorder(biome, -count))  # Order by decreasing count
  
  # Create the stacked bar plot with biome on x-axis
  ggplot(df_filtered, aes(x = biome, fill = level)) +
    geom_bar(position = "stack") +  
    labs(x = "Biome", y = "Number of Regions", fill = "Hierarchy Level") +
    scale_fill_brewer(palette = "Oranges") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

count_biome_levels(regions_with_genus_named)

#


#---------COUNT GCF SPLIT BY GENUS-----------------

#Counts of gcf split by genus
gcf_genus_counts <- plot_bgc_split_counts(bigslice_gcf_id, genus_name, regions_with_genus_named, 1200)
#

#---------COUNT GENUS SPLIT BY GCF-----------------

#Counts of genus split by gcf
genus_gcf_counts <- plot_bgc_split_counts(genus_name, bigslice_gcf_id,regions_with_genus_named, 1000)
#

#---------COUNT TYPE SPLIT BY GENUS-------
type_genus_counts <- plot_bgc_split_counts(type, genus_name, regions_with_genus_named, 100)

#---------COUNT BIOME SPLIT BY GENUS-----------------

#Counts of biome split by genus
biome_genus_counts <- plot_bgc_split_counts(biome, genus_name,regions_with_genus_named, 0)
#

#---------COUNT SUBBIOME SPLIT BY GENUS-----------------

#Counts of subbiome split by genus
subbiome_genus_counts <- plot_bgc_split_counts(subbiome, genus_name,regions_with_genus_named, 0)
#

#---------COUNT INFRABIOME SPLIT BY GENUS-----------------

#Counts of subbiome split by genus
infrabiome_genus_counts <- plot_bgc_split_counts(infrabiome, genus_name,regions_with_genus_named, 0)
#

#---------COUNT MICROBIOME SPLIT BY GENUS-----------------

#Counts of subbiome split by genus
microbiome_genus_counts <- plot_bgc_split_counts(microbiome, genus_name,regions_with_genus_named, 0)
#

#---------COUNT BIOME SPLIT BY GCF---------
#Counts of biome split by gcf
#biome_gcf_counts <- plot_bgc_split_counts(biome, bigslice_gcf_id,regions_with_genus_named, 0, plot=FALSE)
biome_gcf_counts <- plot_bgc_split_counts(biome, bigslice_gcf_id,regions_taxonomy_biome, 0, plot=FALSE)
#

#---------COUNT SUBBIOME SPLIT BY GCF---------
#Counts of biome split by gcf
#subbiome_gcf_counts <- plot_bgc_split_counts(subbiome, bigslice_gcf_id,regions_with_genus_named, 0, plot=FALSE)
subbiome_gcf_counts <- plot_bgc_split_counts(subbiome, bigslice_gcf_id,regions_taxonomy_biome, 0, plot=FALSE)
#

#---------COUNT INFRABIOME SPLIT BY GCF---------
#Counts of biome split by gcf
#infrabiome_gcf_counts <- plot_bgc_split_counts(infrabiome, bigslice_gcf_id,regions_with_genus_named, 0, plot=FALSE)
infrabiome_gcf_counts <- plot_bgc_split_counts(infrabiome, bigslice_gcf_id,regions_taxonomy_biome, 0, plot=FALSE)
#

#---------COUNT MICROBIOME SPLIT BY GCF---------
#Counts of biome split by gcf
#microbiome_gcf_counts <- plot_bgc_split_counts(microbiome, bigslice_gcf_id,regions_with_genus_named, 0, plot=FALSE)
microbiome_gcf_counts <- plot_bgc_split_counts(microbiome, bigslice_gcf_id,regions_taxonomy_biome, 0, plot=FALSE)
#

#---------COUNT BIOME SPLIT BY REGIONS---------
#Counts of biome split by gcf
#biome_region_counts <- plot_bgc_split_counts(biome, region_id,regions_with_genus_named, 0, plot = FALSE)
#

#-----------ASSEMBLY DIVERSITY OF BIOMES-----------

# Create the box plot of Genus Diversity of Assembly by Biome
ggplot(sample_genus_diversity, aes(x = biome, y = gen_assmbly_shannon)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(x = "Biome", y = "Shannon Diversity", title = "Shannon Diversity by Biome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the box plot of Genus Diversity of Assembly by Subbiome
ggplot(sample_genus_diversity, aes(x = subbiome, y = gen_assmbly_shannon)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(x = "Subbiome", y = "Shannon Diversity", title = "Shannon Diversity by Subbiome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



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
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      #theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
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
#

#-------CATEGORY PROCESSING FUNCTIONS-------

# Function to check and relabel
relabel_mixed_categories <- function(x) {
  # Extract all categories
  categories <- unique(unlist(str_extract_all(x, "[A-Za-z]+")))
  
  # Check the number of unique categories
  if (length(categories) > 1) {
    return("{mixed}")
  } else {
    # If only one category, keep the original
    return(x)
  }
}

# Function to clean and reorder a single string
clean_categories <- function(x) {
  #Extract all categories, remove duplicates, and sort
  all_categories <- unique(unlist(str_extract_all(x, "[A-Za-z]+")))
  sorted_categories <- sort(all_categories)
  
  # Combine into a single string with curly braces
  cleaned <- paste0("{", paste(sorted_categories, collapse = ", "), "}")
  return(cleaned)
}

#Analyzing correlation Abundance-Shannon of each Category
calc_category_correlation <- function(df_diversity) {
  # Get unique categories
  unique_categories <- unique(df_diversity$clean_categories)
  
  # Initialize an empty dataframe to store results
  correlation_results <- data.frame(
    category = character(),
    pearson_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each category and calculate Pearson correlation
  for (cat in unique_categories) {
    # Filter data for the current category and remove NA or infinite values
    filtered_data <- df_diversity %>%
      filter(clean_categories == cat) %>%
      filter(is.finite(log10_abundance), is.finite(shannon))
    
    # Check if there are at least 2 data points to calculate correlation
    if (nrow(filtered_data) > 2) {
      # Calculate Pearson correlation
      corr_test <- cor.test(filtered_data$log10_abundance, filtered_data$shannon, method = "pearson")
      
      # Store the results
      correlation_results <- rbind(correlation_results, data.frame(
        category = cat,
        pearson_value = corr_test$estimate,
        p_value = corr_test$p.value
      ))
    }
  }
  
  return((correlation_results))
}


#---------GENUS DIVERSITY FOR GCF-----------

#Diversity values of genus in gcf with richness higher than 1
gcf_genus_diversity <- plot_split_diversity(gcf_genus_counts, bigslice_gcf_id, "genus")

summary(gcf_genus_diversity$r_simpson)
#The distribution is skewed with many low values but a tail extending to high values.
#Some communities are dominated by one species (low Reverse Simpson). Others are more evenly distributed (high Reverse Simpson).
summary(gcf_genus_diversity$shannon)
#The distribution shows a wide range of diversity: from monocultures (0.0000) to highly diverse communities (5.4264).
#The relatively higher values in the Shannon Index compared to Reverse Simpson suggest that rare species are contributing to diversity in some communities.

#Merging abundance information to diversity
gcf_genus_diversity$log10_abundance <- log10(rowSums(gcf_genus_counts[,-1]))

#Combine product categories and create new columns
combined_categories <- regions_with_genus_named %>%
  group_by(bigslice_gcf_id) %>%
  summarise(
    adjusted_product_categories = paste0(
      sort(unique(unlist(str_extract_all(adjusted_product_categories, "\\{[^}]+\\}")))), 
      collapse = ", "), 
    .groups = "drop"
  ) %>%
  mutate(
    clean_categories = sapply(adjusted_product_categories, clean_categories),
    relabel_mixed_categories = sapply(adjusted_product_categories, relabel_mixed_categories)
  )

# Add the new columns to gcf_genus_diversityy
gcf_genus_diversity <- gcf_genus_diversity %>%
  left_join(combined_categories %>% select(bigslice_gcf_id, clean_categories, relabel_mixed_categories),
            by = "bigslice_gcf_id")

#Clean categories
print(length(unique(gcf_genus_diversity$clean_categories)))

#Product category pie char
pie(table(gcf_genus_diversity$clean_categories)[table(gcf_genus_diversity$clean_categories) > 200], main = "Pie Chart of Clean Categories (> 200 counts)")
# Increase margin size
#par(mar=c(14,4,4,4))
#barplot(sort(table(gcf_genus_diversity$clean_categories), decreasing = TRUE), las=2, col=rgb(0.2,0.4,0.6,0.6), main="Category distribution", ylab ="Number of regions")
plot_bgc_split_counts(clean_categories, bigslice_gcf_id, gcf_genus_diversity, 0, split=FALSE, color = "pink2")

#Scatterplot representing influence of abundance on shannon index
ggplot(gcf_genus_diversity, aes(x=log10_abundance, y=shannon)) + 
  geom_point(shape = 19, alpha=0.5, color = "coral") 

ggplot(gcf_genus_diversity[gcf_genus_diversity$clean_categories %in% c("{NRPS}", "{PKS}","{RiPP}","{terpene}"), ], aes(x=log10_abundance, y=shannon , color = clean_categories)) + 
  geom_point(shape = 19, alpha=0.7) 

ggplot(gcf_genus_diversity[gcf_genus_diversity$clean_categories %in% c("{NRPS}", "{PKS}", "{RiPP}", "{terpene}"), ], 
       aes(x = log10_abundance, y = shannon, color = clean_categories)) + 
  geom_point(shape = 19, alpha = 1) +
  facet_wrap(~ clean_categories, scales = "free") +
  labs(x = "Log10 Abundance", y = "Shannon Index", color = "Category") +
  theme_minimal() +
  theme(legend.position = "none")

#Calculate correlation of all categories with log10_abundance
full_cat_correlation <- calc_category_correlation(gcf_genus_diversity)

ggplot(data=full_cat_correlation[full_cat_correlation$p_value < .05,], aes(x = reorder(category, pearson_value), y = pearson_value, fill = "coral")) +
  geom_col(position = position_stack(), ) +
  geom_hline(yintercept = .2, color = "aquamarine3") +
  geom_hline(yintercept = .4, color = "aquamarine2") +
  geom_hline(yintercept = .6, color = "aquamarine1") +
  ylab("Pearson Correlation Index") + 
  xlab("Category") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = -75, hjust = 0, vjust = 1))

#---------TYPE OF TOP GCF-GENUS DIVERSE-------


#Isolating gcfs with high top diversity values
top_gcf_genus_diversity <- gcf_genus_diversity[gcf_genus_diversity$shannon >= 3,] # = 495 of 18383
#top_gcf_genus_diversity <- gcf_genus_diversity[gcf_genus_diversity$shannon >= 1,] # = 6303 of 18383
top_gcf_genus_diversity <- top_gcf_genus_diversity %>%
  left_join(regions_with_genus_named %>% select(bigslice_gcf_id, type),
            by = "bigslice_gcf_id")
top_gcf_genus_diversity <- top_gcf_genus_diversity[!duplicated(top_gcf_genus_diversity), ]
top_gcf_genus_diversity <- top_gcf_genus_diversity[!top_gcf_genus_diversity$bigslice_gcf_id == -1, ]
top_gcf_genus_diversity <- top_gcf_genus_diversity %>%
  arrange(desc(shannon)) #%>% # arrange in descending order 
  #slice(1:100)

# Step 1: Concatenate all type values per bigslice_gcf_id and shannon
summ_top_gcf_genus_diversity <- top_gcf_genus_diversity %>%
  group_by(bigslice_gcf_id, shannon) %>%
  summarise(type = paste(type, collapse = ","), .groups = "drop")

# Step 2: Remove repeated words
summ_top_gcf_genus_diversity <- summ_top_gcf_genus_diversity %>%
  mutate(type = sapply(strsplit(type, ","), function(x) paste(unique(x), collapse = ","))) %>%
  arrange(desc(shannon))

ggplot(summ_top_gcf_genus_diversity[c(1:14, 16:18),], aes(x = reorder(bigslice_gcf_id, -shannon))) +
  geom_point(aes(y =  shannon, color = type), size = 3) + 
  labs(x = "Type", y = "Shannon Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))

#

#---------GENUS DIVERSITY OF GCF-------

#Diversity values of genus in gcf with richness higher than 1
gcf_genus_diversity <- plot_split_diversity(gcf_genus_counts, genus_name, "gcf")

summary(gcf_genus_diversity$r_simpson)
#On average very low
summary(gcf_genus_diversity$shannon)
#On average very low

boxplot(gcf_genus_diversity$richness, ylim = c(0,200))
barplot(sort(gcf_genus_diversity$richness, decreasing = TRUE))

#Calculating gcf_genus-diversity by biome
gcf_genus_diversity_by_biome <- multiple_diversity(regions_with_genus_named, "biome")

#Calculating gcf_genus-diversity by subbiome
gcf_genus_diversity_by_subbiome <- multiple_diversity(regions_with_genus_named, "subbiome")

#Calculating gcf_genus-diversity by subbiome
gcf_genus_diversity_by_infrabiome <- multiple_diversity(regions_with_genus_named, "infrabiome")

#Calculating gcf_genus-diversity by subbiome
gcf_genus_diversity_by_microbiome <- multiple_diversity(regions_with_genus_named, "microbiome")

#---------GENUS DIVERSITY FOR BIOME---------

#Diversity values of genus in biome
biome_genus_diversity <- plot_split_diversity(biome_genus_counts, biome, "genus")

summary(biome_genus_diversity$shannon)
#

#---------GENUS DIVERSITY FOR SUBBIOME---------

#Diversity values of genus in biome
subbiome_genus_diversity <- plot_split_diversity(subbiome_genus_counts, subbiome, "genus")

summary(subbiome_genus_diversity$shannon)
#

#---------GENUS DIVERSITY FOR INFRABIOME---------

#Diversity values of genus in biome
infrabiome_genus_diversity <- plot_split_diversity(infrabiome_genus_counts, infrabiome, "genus")

summary(infrabiome_genus_diversity$shannon)
#

#---------GENUS DIVERSITY FOR MICROBIOME---------

#Diversity values of genus in biome
microbiome_genus_diversity <- plot_split_diversity(microbiome_genus_counts, microbiome, "genus")

summary(infrabiome_genus_diversity$shannon)
#

#-----------COMPARING GENUS DIVERSITY FOR BIOME LEVELS--------

# Create a combined dataframe
biome_level_genus_richness <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_genus_diversity$richness),
  data.frame(Level = "Subbiome", Shannon = subbiome_genus_diversity$richness),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_genus_diversity$richness),
  data.frame(Level = "Microbiome", Shannon = microbiome_genus_diversity$richness)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

# Create the boxplot
ggplot(biome_level_genus_richness, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Genus Richness") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_genus_diversity$richness, biome_genus_diversity$richness)
# Kolmogorov-Smirnov Test: D = 0.31515, p-value = 0.2082
wilcox.test(biome_genus_diversity$richness, biome_genus_diversity$richness)
# Wilcoxon Rank-Sum Test: W = 240, p-value = 0.8763
# Test indicate distributions are not significantly different

# Create a combined dataframe
biome_level_genus_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_genus_diversity$shannon),
  data.frame(Level = "Subbiome", Shannon = subbiome_genus_diversity$shannon),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_genus_diversity$shannon),
  data.frame(Level = "Microbiome", Shannon = microbiome_genus_diversity$shannon)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

#boxplot(biome_genus_diversity$shannon, subbiome_genus_diversity$shannon, infrabiome_genus_diversity$shannon, microbiome_genus_diversity$shannon, xlab="Genus Shannon Index for each Biome Level")

# Create the boxplot
ggplot(biome_level_genus_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greens") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "Genus Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

#Analyzing if distributions are significantly different
ks.test(biome_genus_diversity$shannon,microbiome_genus_diversity$shannon)
# Kolmogorov-Smirnov Test: D = 0.25455, p-value = 0.4388
wilcox.test(biome_genus_diversity$shannon,microbiome_genus_diversity$shannon)
# Wilcoxon Rank-Sum Test: W = 221, p-value = 0.5671
# Test indicate distributions are not significantly different

#-------SHANNON DIVERSITY COMPARISON FUNCTION-----

#Function to compare the shannon diversity assembly-genus with genus-gcf per bioem
plot_triple_comparison <- function(sample_genus_diversity, gcf_genus_diversity, level, plot = 2) {
  
  # Capture the input as symbols for tidy evaluation
  level_sym <- ensym(level)
  
  # Calculate the mean of gen_assmbly_shannon for each biome to determine the order
  biome_order <- sample_genus_diversity %>%
    group_by(!!level_sym) %>%
    summarize(mean_shannon = mean(gen_assmbly_shannon, na.rm = TRUE)) %>%
    arrange(mean_shannon) %>%
    pull(!!level_sym)
  
  # Combine the two dataframes into one for plotting
  combined_df <- bind_rows(
    sample_genus_diversity %>% 
      select(shannon = gen_assmbly_shannon, !!level_sym) %>%
      mutate(source = "genus_assembly"),
    gcf_genus_diversity %>% 
      select(shannon, !!level_sym) %>%
      mutate(source = "gcf_genus")
  )
  
  if (plot == 1) {
    
    # Plot the boxplots
    ggplot(combined_df, aes(x = factor(!!level_sym, levels = biome_order), y = shannon, fill = source)) +
      geom_boxplot() +
      labs( x = deparse1(substitute(level)), y = "Shannon Diversity", title = paste("Shannon Diversity by", deparse1(substitute(level)) ), fill = "Source") +
      scale_fill_manual(values = c("genus_assembly" = "hotpink1", "gcf_genus" = "darkslategray3")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = -45, hjust = 0))
    
  } else {
    
    # Calculate means and standard errors
    mean_se_df <- combined_df %>%
      group_by(!!level_sym, source) %>%
      summarize(
        mean_shannon = mean(shannon, na.rm = TRUE),
        se = sd(shannon, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    # Plot the means with standard error bars and lines
    ggplot(mean_se_df, aes(x = factor(!!level_sym, levels = biome_order), y = mean_shannon, color = source, group = source)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_shannon - se, ymax = mean_shannon + se), width = 0.2) +
      geom_line() +
      labs(
        x = deparse1(substitute(level)),
        y = "Shannon Diversity (Mean ± SE)",
        title = paste("Mean Shannon Diversity by", deparse1(substitute(level)), "with Standard Error"),
        color = "Source"
      ) +
      scale_color_manual(values = c("genus_assembly" = "hotpink1", "gcf_genus" = "darkslategray3")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = -45, hjust = 0))
  }
}

# Define the function to calculate Spearman correlation of means
calculate_spearman_correlation <- function(sample_df, gcf_df, level, meth = "spearman") {
  
  # Capture the input as symbols for tidy evaluation
  level_sym <- ensym(level)
  
  # Calculate the mean of gen_assmbly_shannon for each biome
  sample_means <- sample_df %>%
    group_by(!!level_sym) %>%
    summarize(mean_gen_assmbly = mean(gen_assmbly_shannon, na.rm = TRUE), .groups = "drop")
  
  # Calculate the mean of shannon for each biome
  gcf_means <- gcf_df %>%
    group_by(!!level_sym) %>%
    summarize(mean_gcf_genus = mean(shannon, na.rm = TRUE), .groups = "drop")
  
  # Merge the means by biome
  combined_means <- inner_join(sample_means, gcf_means, by = deparse1(substitute(level)))
  
  # Calculate the Spearman correlation
  spearman_corr <- cor.test(combined_means$mean_gen_assmbly, combined_means$mean_gcf_genus, method = meth)
  
  print(spearman_corr)
}

#-----------COMPARING GENUS-GCF DIVERSITY WITH BIOME-GENUS DIVERSITY----------

#Boxplot
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_biome,  biome, plot =1)
#Mean error
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_biome,  biome)

# Calculating Correlation
calculate_spearman_correlation(sample_genus_diversity, gcf_genus_diversity_by_biome,  biome)
#S = 3596, p-value = 0.05963, rho -0.3830769 
#Moderate negative correlation, but it is not statistically significant.

#-----------COMPARING GENUS-GCF DIVERSITY WITH SUBBIOME-GENUS DIVERSITY----------

#Boxplot
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_subbiome,  subbiome, plot =1)
#Mean error
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_subbiome,  subbiome)

# Calculating Spearman Correlation
calculate_spearman_correlation(sample_genus_diversity, gcf_genus_diversity_by_subbiome,  subbiome)
#S = 12802, p-value = 0.06789, rho -0.295749 
#Weak negative correlation, but it is not statistically significant.

#-----------COMPARING GENUS-GCF DIVERSITY WITH INFRABIOME-GENUS DIVERSITY----------

#Boxplot
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_infrabiome,  infrabiome, plot =1)
#Mean error
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_infrabiome,  infrabiome)

# Calculating Correlation
calculate_spearman_correlation(sample_genus_diversity, gcf_genus_diversity_by_infrabiome,  infrabiome)
#S = 31100, p-value = 0.003239, rho = -0.4072398
#Strong negative correlation and statistically significant.

#-----------COMPARING GENUS-GCF DIVERSITY WITH MICROBIOME-GENUS DIVERSITY----------

#Boxplot
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_microbiome,  microbiome, plot =1)
#Mean error
plot_triple_comparison(sample_genus_diversity, gcf_genus_diversity_by_microbiome,  microbiome)

# Calculating Spearman Correlation
calculate_spearman_correlation(sample_genus_diversity, gcf_genus_diversity_by_microbiome,  microbiome)
#S = 10910, p-value = 0.001307, rho -0.5280112 
#Strong negative correlation and statistically significant.
calculate_spearman_correlation(sample_genus_diversity, gcf_genus_diversity_by_microbiome,  microbiome, meth = "pearson")

#---------GCF DIVERSITY FOR BIOME-------
#GCF Diversity values for biome
biome_gcf_diversity <- plot_split_diversity(biome_gcf_counts, biome, "genus")
#

#---------GCF DIVERSITY FOR SUBBIOME-------
#GCF Diversity values for subbiome
subbiome_gcf_diversity <- plot_split_diversity(subbiome_gcf_counts, subbiome, "genus")
#

#---------GCF DIVERSITY FOR INFRABIOME-------
#GCF Diversity values for infrabiome
infrabiome_gcf_diversity <- plot_split_diversity(infrabiome_gcf_counts, infrabiome, "genus")
#

#---------GCF DIVERSITY FOR MICROBIOME-------
#GCF Diversity values for microbiome
microbiome_gcf_diversity <- plot_split_diversity(microbiome_gcf_counts, microbiome, "genus")
#

#-----------COMPARING GCF DIVERSITY FOR BIOME LEVELS  --------
# Create a combined dataframe
biome_level_gcf_shannon <- bind_rows(
  data.frame(Level = "Biome", Shannon = biome_gcf_diversity$shannon),
  data.frame(Level = "Subbiome", Shannon = subbiome_gcf_diversity$shannon),
  data.frame(Level = "Infrabiome", Shannon = infrabiome_gcf_diversity$shannon),
  data.frame(Level = "Microbiome", Shannon = microbiome_gcf_diversity$shannon)
)%>%
  mutate(Level = factor(Level, levels = c("Biome", "Subbiome", "Infrabiome", "Microbiome")))  # Define order

#boxplot(biome_gcf_diversity$shannon, subbiome_gcf_diversity$shannon, infrabiome_gcf_diversity$shannon, microbiome_gcf_diversity$shannon, xlab="GCF Shannon Index for each Biome Level")

# Create the boxplot
ggplot(biome_level_gcf_shannon, aes(x = Level, y = Shannon, fill = Level)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Reds") +  # Try "Set1", "Set2", "Pastel1", "Dark2", etc
  labs(x = "Biome Level", y = "GCF Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

ks.test(biome_gcf_diversity$shannon,microbiome_gcf_diversity$shannon)
# Kolmogorov-Smirnov Test: D = 0.27879, p-value = 0.3306
#D = 0.31439, p-value = 0.1959
wilcox.test(biome_gcf_diversity$shannon,microbiome_gcf_diversity$shannon)
# Wilcoxon Rank-Sum Test: W = 254, p-value = 0.895
# W = 286, p-value = 0.65
# Test indicate distributions are not significantly different

#-------COMPARING FUNCTION FOR BIOLEVELS-GCF WITH BIOLEVELS-GENUS DIVERSITY--------------
plot_biolevel_gcf_genus <- function(level_gcf_diversity, level_genus_diversity, level) {
  
  # Merge the dataframes on the selected 'level' column
  merged_data <- level_gcf_diversity %>%
    inner_join(level_genus_diversity, by = rlang::as_name(ensym(level)), suffix = c("_gcf", "_genus")) %>%
    arrange(richness_genus)  # Sorting by increasing richness of biome_genus
  
  
  # Create the plot
  ggplot(merged_data, aes(x = reorder(!!ensym(level), richness_genus))) +
    geom_point(aes(y = richness_gcf, color = "GCF"), size = 3) + 
    geom_point(aes(y = richness_genus, color = "Genus"), size = 3) + 
    geom_line(aes(y = richness_gcf, group = 1, color = "GCF"), linetype = "dashed") +
    geom_line(aes(y = richness_genus, group = 1, color = "Genus"), linetype = "solid") +
    #geom_errorbar(aes(ymin = richness_gcf - sd(richness_gcf), ymax = richness_gcf + sd(richness_gcf), color = "GCF"), width = 0.2) +
    #geom_errorbar(aes(ymin = richness_genus - sd(richness_genus), ymax = richness_genus + sd(richness_genus), color = "Genus"), width = 0.2) +
    labs(x = as_label(ensym(level)), y = "Richness", color = "Diversity Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -60, hjust = 0))
}
#

#-----------COMPARING BIOME-GCF RICHNESS WITH BIOME-GENUS DIVERSITY--------

plot_biolevel_gcf_genus(biome_gcf_diversity, biome_genus_diversity, biome)

cor.test(biome_gcf_diversity[-16,]$richness, biome_genus_diversity$richness, method = "spearman")

#

#-----------COMPARING SUBBIOME-GCF RICHNESS WITH SUBBIOME-GENUS DIVERSITY--------

plot_biolevel_gcf_genus(subbiome_gcf_diversity[!(is.na(subbiome_gcf_diversity$subbiome)), ], 
                        subbiome_genus_diversity[!(is.na(subbiome_genus_diversity$subbiome)), ], subbiome)

cor.test(subbiome_gcf_diversity$richness, subbiome_genus_diversity$richness, method = "spearman")
#rho: .933

#Splitting data between enviromental and host associated
env_subbiomes <- c( "Aquatic:Aquaculture" , "Aquatic:Estuary" , "Aquatic:Freshwater", "Aquatic:Lentic" , "Aquatic:Marine" , 
                    "Aquatic:Non-marine Saline and Alkaline", "Aquatic:Sediment" , "Aquatic:Thermal springs","Terrestrial:Soil" )
host_subbiomes <- c("Algae:Red algae" , "Animal:Digestive system" , "Arthropoda:Digestive system", "Birds:Digestive system" , 
                    "Fish:Digestive system" , "Human:Digestive system"  , "Human:Skin" , "Insecta:Digestive system" , 
                    "Mammals:Digestive system" , "Mammals:Gastrointestinal tract" , "Mammals:Respiratory system" , 
                    "Microbial:Bacteria" , "Plants:Phylloplane" , "Plants:Rhizosphere" ,"Plants:Root"  )

plot_biolevel_gcf_genus(subbiome_gcf_diversity[subbiome_gcf_diversity$subbiome %in% env_subbiomes, ], 
                        subbiome_genus_diversity[subbiome_genus_diversity$subbiome %in% env_subbiomes, ], subbiome)
cor.test(subbiome_gcf_diversity[subbiome_gcf_diversity$subbiome %in% env_subbiomes, ]$richness, subbiome_genus_diversity[subbiome_genus_diversity$subbiome %in% env_subbiomes, ]$richness, method = "spearman")
#rho: 0.9333


plot_biolevel_gcf_genus(subbiome_gcf_diversity[subbiome_gcf_diversity$subbiome %in% host_subbiomes, ], 
                        subbiome_genus_diversity[subbiome_genus_diversity$subbiome %in% host_subbiomes, ], subbiome)
cor.test(subbiome_gcf_diversity[subbiome_gcf_diversity$subbiome %in% host_subbiomes, ]$richness, subbiome_genus_diversity[subbiome_genus_diversity$subbiome %in% host_subbiomes, ]$richness, method = "spearman")
#rho: 0.9722

#

#-----------COMPARING INFRABIOME-GCF RICHNESS WITH INFRABIOME-GENUS DIVERSITY--------

plot_biolevel_gcf_genus(infrabiome_gcf_diversity[!(is.na(infrabiome_gcf_diversity$infrabiome)), ], 
                        infrabiome_genus_diversity[!(is.na(infrabiome_genus_diversity$infrabiome)), ], infrabiome)

cor.test(infrabiome_gcf_diversity[-47,]$richness, infrabiome_genus_diversity$richness, method = "spearman")
#

#-----------COMPARING MICROBIOME-GCF RICHNESS WITH MICROBIOME-GENUS DIVERSITY--------

plot_biolevel_gcf_genus(microbiome_gcf_diversity[!(is.na(microbiome_gcf_diversity$microbiome)), ],
                        microbiome_genus_diversity[!(is.na(microbiome_gcf_diversity$microbiome)), ], microbiome)

cor.test(microbiome_gcf_diversity$richness, microbiome_genus_diversity$richness, method = "spearman")
#

#-------INDICATOR FUNCTIONS-------------

#Calculates abundance of var in biome
calc_abundance <- function(df,biome_var, var) {
  # Count occurrences of each GCF per assembly
  abundance_df <- df %>%
    count(assembly, !!sym(var)) %>%  # Dynamically use gcf_var
    pivot_wider(names_from = !!sym(var), values_from = n, values_fill = list(n = 0))
  
  # Ensure only unique assembly-biome pairs are kept
  df_unique <- df %>%
    select(assembly, !!sym(biome_var)) %>%
    distinct()
  
  # Merge the biome variable into the abundance dataframe
  out_df <- df_unique %>%
    left_join(abundance_df, by = "assembly")
  
  return(out_df)
}

#Calculates total specificty. Takes abundancae dataframe as input, grouping biome and column to compare.
calc_total_specificity <- function(abundance_df, grouping_var, var) {
  
  # Convert abundance data to long format
  abundance_long <- abundance_df %>%
    pivot_longer(cols = -all_of(c("assembly", grouping_var)), 
                 names_to = var, values_to = "abundance")
  
  # Calculate total abundance per genus per grouping variable
  total_per_group <- abundance_long %>%
    group_by(across(all_of(c(grouping_var, var)))) %>%
    summarise(total_abundance_group = sum(abundance, na.rm = TRUE), .groups = "drop")
  
  # Calculate total abundance per genus across all groups
  total_overall <- abundance_long %>%
    group_by(across(all_of(var))) %>%
    summarise(total_abundance_overall = sum(abundance, na.rm = TRUE), .groups = "drop")
  
  # Merge and compute specificity
  specificity_df <- total_per_group %>%
    left_join(total_overall, by = var) %>%
    mutate(specificity = total_abundance_group / total_abundance_overall)
  
  # Reshape to wide format: biomes as rows, genera as columns
  specificity_wide <- specificity_df %>%
    select(!!sym(grouping_var), !!sym(var), specificity) %>%  # Keep only relevant columns
    pivot_wider(names_from = !!sym(var), values_from = specificity, values_fill = list(specificity = 0))  # Reshape
  
  return(specificity_df)
}

#Reshape specificity
reshape_table <- function(df, biome_var, var, value_var) {
  df %>%
    select(!!sym(biome_var), !!sym(var), !!sym(value_var)) %>%  # Keep only relevant columns
    pivot_wider(names_from = !!sym(var), values_from = !!sym(value_var), 
                values_fill = setNames(list(0), value_var))  # Corrected values_fill syntax
}

#Function for calculating fidelity
calc_fidelity <- function(abundance_df, biome_var, var) {
  # Convert to long format
  abundance_long <- abundance_df %>%
    pivot_longer(cols = -c(assembly, !!sym(biome_var)), names_to = var, values_to = "abundance")
  
  # Count assemblies where abundance > 0 per genus in each biome
  presence_counts <- abundance_long %>%
    group_by(!!sym(biome_var), !!sym(var)) %>%
    summarise(present_count = sum(abundance > 0, na.rm = TRUE), .groups = "drop")
  
  # Count total number of assemblies per biome
  total_assemblies <- abundance_long %>%
    select(assembly, !!sym(biome_var)) %>%
    distinct() %>%
    group_by(!!sym(biome_var)) %>%
    summarise(total_assemblies = n(), .groups = "drop")
  
  # Compute fidelity (present_count / total assemblies)
  fidelity_df <- presence_counts %>%
    left_join(total_assemblies, by = biome_var) %>%
    mutate(fidelity = present_count / total_assemblies) %>%
    select(!!sym(biome_var), !!sym(var), fidelity)  # Keep only relevant columns
  
  # Reshape to wide format: biomes as rows, genera as columns
  fidelity_wide <- fidelity_df %>%
    pivot_wider(names_from = !!sym(var), values_from = fidelity, values_fill = 0)
  
  return(fidelity_df)
  #return(fidelity_wide)
}

#Pipeline to get indicator table
get_indicator <- function(df, biome_var, var, abundance=FALSE) {
  print("Step 1: Calculate Genus Abundance")
  if (abundance){
    abundance <- df
  }else{
    abundance <- calc_abundance(df, biome_var, var)
  }
  
  
  print("Step 2: Calculate Total Specificity for Genus on Biome")
  total_specificity <- calc_total_specificity(abundance, biome_var, var)
  
  # Step 3: Reshape Specificity Table
  #genus_specificity <- reshape_specificity_table(genus_total_specificity, biome_var, var, column_var)
  
  print("Step 4: Calculate Genus Fidelity")
  fidelity <- calc_fidelity(abundance, biome_var, var)
  
  print("Step 5: Perform Indicator Species Analysis (ISA)")
  indicator <- cbind(fidelity, total_specificity[,5])
  indicator$indicator <- indicator$fidelity*indicator$specificity*100
  
  
  # Return all results as a list
  return(indicator)
  #return(list(genus_abundance = genus_abundance,genus_specificity = genus_specificity,genus_fidelity = genus_fidelity,genus_indicator = genus_indicator))
}


#---------GCF SPECIFICITY FOR GENUS---------

#How specific the gcf is to the genus
genus_gcf_specificity <- read.delim(file = "atlas_tables/genus_gcf_specificity.csv", sep = ",", header=FALSE)
colnames(genus_gcf_specificity) <- c("genus_name","gfc_id","abundance","specificity")

#Picking highest GCF specificity for a GENUS
genus_gcf_top_specificity <- genus_gcf_specificity %>%
  group_by(gfc_id) %>%           # Group by gfc_id
  slice_max(order_by = specificity, n = 1) %>%  # Select row with the highest specificity
  ungroup()                      # Ungroup the data
#Summary statistics of highest specificity of a gcf to a genus combination
summary(genus_gcf_top_specificity$specificity)

#par(mfrow =c(1,1))
#boxplot(genus_gcf_top_specificity$specificity, xlab="Top Genus Specificity for each GCF")

# Create the boxplot
ggplot(genus_gcf_top_specificity, aes( y = specificity)) +
  geom_boxplot(fill="slateblue", alpha=0.4, width=0.7) +
  xlim(-1,1) +
  labs( y = "Top GCF Specificity") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())


extract_gcf <- function(top_specificity, counts, filter_condition) {
  # Apply the filter condition dynamically
  filtered_gcf <- subset(top_specificity, eval(parse(text = filter_condition)))
  
  # Extract unique GCF IDs from the filtered data
  outlier_gcf <- unique(filtered_gcf$gfc_id)
  cat(length(outlier_gcf), "GCFs meet the filter condition\n")
  
  # Calculate total counts for each GCF in outliers
  outlier_counts <- counts[counts$bigslice_gcf_id %in% outlier_gcf, ]
  outlier_total_counts <- rowSums(outlier_counts[,-1])
  
  # Return total counts
  return(outlier_total_counts)
}

#Extract and analyze gcf below the 1st quartile. (specificity lower than 0.4)
outlier_gcf_total_counts <- extract_gcf(genus_gcf_top_specificity, genus_gcf_counts, "specificity < 0.4")
summary(outlier_gcf_total_counts)

#Extract and analyze gcf above the 1st quartile. (specificity higher than 0.4)
main_gcf_total_counts <- extract_gcf(genus_gcf_top_specificity, genus_gcf_counts, "specificity > 0.4")
summary(main_gcf_total_counts)

#Analyzing if distributions are significantly different
ks.test(main_gcf_total_counts, outlier_gcf_total_counts)
# Kolmogorov-Smirnov Test: D = 0.34409, p-value < 2.2e-16
# The very low p-value indicates that the distributions are significantly different.
wilcox.test(main_gcf_total_counts, outlier_gcf_total_counts)
# Wilcoxon Rank-Sum Test: W = 13234857, p-value < 2.2e-16
# This also shows a significant difference in the central tendencies of the two distributions.

# Combine vectors into a data frame
top_gcf_distributions <- data.frame(
  value = c(main_gcf_total_counts, outlier_gcf_total_counts),
  group = c(rep("Above 1st Qu.", length(main_gcf_total_counts)), 
            rep("Below 1st Qu.", length(outlier_gcf_total_counts)))
  )
top_gcf_distributions$log_value <- log10(top_gcf_distributions$value + 1)

ks.test(main_gcf_total_counts, outlier_gcf_total_counts)
# Kolmogorov-Smirnov Test: D = 0.34409, p-value < 2.2e-16
# The very low p-value indicates that the distributions are significantly different.
wilcox.test(main_gcf_total_counts, outlier_gcf_total_counts)
# Wilcoxon Rank-Sum Test: W = 13234857, p-value < 2.2e-16
# This also shows a significant difference in the central tendencies of the two distributions.


# Violin plot with log-transformed data
ggplot(top_gcf_distributions, aes(x = group, y = log_value, fill = group)) + 
  geom_violin(trim = FALSE, alpha = 0.5) + 
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) + 
  labs(
    title = "Violin Plot of Log-Transformed GCF Total Counts",
    x = "Top GCF Specificity to Genus",
    y = "Log10( GCF Total Counts + 1 )"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#GCFs with lower specificity values are on average more widespread


#Scatterplot representing influence of abundance on specificity
ggplot(genus_gcf_specificity, aes(y=specificity, x=log10(abundance))) + 
  geom_point(shape = 19, alpha=0.5, color = "coral") 

#Scatterplot representing influence of abundance on top specificity
ggplot(genus_gcf_top_specificity, aes(y=specificity, x=log10(abundance))) + 
  geom_point(shape = 19, alpha=0.5, color = "coral") 
