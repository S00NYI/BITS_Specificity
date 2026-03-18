################################################################################
## eCLIP Cellular Specificity and Variation Sensitivity Analysis:
## Written by Soon Yi with Antigravity
## Created: 2026-03-17
## Last Edited: 2026-03-17
## Figure 6
################################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(ggsignif)

## Set up basic parameters:
################################################################################
baseDir = '~/Repos/BITS_Specificity/Dataset/Analysis/'
K = 5
DATE = '20260316'

# Load eCLIP IS/VS summary (from eCLIP_motif_analysis.R):
eCLIP_summary = read_csv(paste0(baseDir, 'eCLIP_all/output/', DATE, '/eCLIP_', K, 'mer_IS_VS_summary.csv'),
                         col_names = TRUE, show_col_types = FALSE)
eCLIP_summary = data.frame(eCLIP_summary)

# Load RBNS IS/VS metrics:
RBNS_metrics = read_csv(paste0(baseDir, 'RBNS/', K, 'mer_metrics.csv'),
                        col_names = TRUE, show_col_types = FALSE)
RBNS_metrics = data.frame(RBNS_metrics)
colnames(RBNS_metrics) = c('RBP', 'IS', 'VS')

# Define unconventional RBPs:
unconventional_RBPs = c('AATF', 'APOBEC3C', 'BCCIP', 'CDC40', 'EIF3D', 'EIF3H',
                        'FKBP4', 'GNL3', 'GRWD1', 'GTF2F1', 'MTPAP', 'NIP7',
                        'NPM1', 'PHF6', 'PUS1', 'RPS11', 'RPS3', 'SBDS',
                        'SF3A3', 'SLBP', 'SMNDC1', 'SUB1', 'UCHL5', 'UTP3',
                        'YWHAG', 'ZNF622')
################################################################################

## Compute representative CS and C-VS for each RBP:
################################################################################
# For RBPs with both cell lines, use averaged value.
# For RBPs with only one cell line, use that value.
eCLIP_summary$CS = ifelse(!is.na(eCLIP_summary$IS_avg), eCLIP_summary$IS_avg,
                   ifelse(!is.na(eCLIP_summary$IS_K562), eCLIP_summary$IS_K562,
                          eCLIP_summary$IS_HepG2))

eCLIP_summary$CVS = ifelse(!is.na(eCLIP_summary$VS_avg), eCLIP_summary$VS_avg,
                    ifelse(!is.na(eCLIP_summary$VS_K562), eCLIP_summary$VS_K562,
                           eCLIP_summary$VS_HepG2))
################################################################################

## Panel 1A: K562 vs HepG2 Cellular Specificity:
################################################################################
# Filter for RBPs with data in both cell lines:
both_cells = eCLIP_summary %>% filter(!is.na(IS_K562) & !is.na(IS_HepG2))

Corr_CS = cor(both_cells$IS_K562, both_cells$IS_HepG2, method = 'pearson')

ggplot(both_cells, aes(x = IS_K562, y = IS_HepG2)) +
  geom_point(fill = 'cornflowerblue', alpha = 1, shape = 21, size = 4) +
  geom_smooth(method = 'lm', linetype = 'solid', color = 'cornflowerblue',
              se = FALSE, linewidth = 0.5, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
  labs(x = 'K562 Cellular Specificity',
       y = 'HepG2 Cellular Specificity',
       title = 'K562 vs HepG2 Cellular Specificity') +
  scale_x_continuous(trans = 'log2', limits = c(1, 33)) +
  scale_y_continuous(trans = 'log2', limits = c(1, 33)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  annotate('text', x = Inf, y = 1, label = sprintf('R = %.2f', Corr_CS),
           hjust = 1.1, vjust = 1.5, size = 5, color = 'cornflowerblue')
################################################################################

## Panel 1B: K562 vs HepG2 Cellular Variation Sensitivity:
################################################################################
Corr_CVS = cor(both_cells$VS_K562, both_cells$VS_HepG2, method = 'pearson')

ggplot(both_cells, aes(x = VS_K562, y = VS_HepG2)) +
  geom_point(fill = 'darkgoldenrod1', alpha = 1, shape = 21, size = 4) +
  geom_smooth(method = 'lm', linetype = 'solid', color = 'darkgoldenrod1',
              se = FALSE, linewidth = 0.5, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50') +
  labs(x = 'K562 Cellular Variation Sensitivity',
       y = 'HepG2 Cellular Variation Sensitivity',
       title = 'K562 vs HepG2 Cellular Variation Sensitivity') +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  annotate('text', x = Inf, y = -Inf, label = sprintf('R = %.2f', Corr_CVS),
           hjust = 1.1, vjust = -0.5, size = 5, color = 'darkgoldenrod1')
################################################################################

## Panel 2: Fold Change CS/IS and C-VS/in vitro VS:
################################################################################
# Merge eCLIP and RBNS data for RBPs present in both:
eCLIP_RBNS = merge(eCLIP_summary[, c('RBP', 'CS', 'CVS')],
                   RBNS_metrics,
                   by = 'RBP',
                   suffixes = c('_eCLIP', '_RBNS'))

# Compute fold changes:
eCLIP_RBNS$FC_CS_IS = eCLIP_RBNS$CS / eCLIP_RBNS$IS
eCLIP_RBNS$FC_CVS_VS = eCLIP_RBNS$CVS / eCLIP_RBNS$VS

# Order by RBNS IS (high to low):
eCLIP_RBNS = eCLIP_RBNS %>% arrange(desc(IS))
eCLIP_RBNS$RBP = factor(eCLIP_RBNS$RBP, levels = eCLIP_RBNS$RBP)

## Panel 2A: CS / IS Fold Change:
################################################################################
ggplot(eCLIP_RBNS, aes(x = RBP, y = FC_CS_IS)) +
  geom_point(fill = 'cornflowerblue', alpha = 1, shape = 21, size = 4) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50') +
  labs(x = 'RBP (ordered by Inherent Specificity)',
       y = 'Fold Change (CS / IS)',
       title = 'Cellular vs Inherent Specificity') +
  scale_y_continuous(trans = 'log2') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
################################################################################

## Panel 2B: C-VS / in vitro VS Fold Change:
################################################################################
eCLIP_RBNS_VS = eCLIP_RBNS %>% arrange(desc(VS))
eCLIP_RBNS_VS$RBP = factor(eCLIP_RBNS_VS$RBP, levels = eCLIP_RBNS_VS$RBP)

ggplot(eCLIP_RBNS_VS, aes(x = RBP, y = FC_CVS_VS)) +
  geom_point(fill = 'darkgoldenrod1', alpha = 1, shape = 21, size = 4) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50') +
  labs(x = 'RBP (ordered by in vitro VS)',
       y = 'Fold Change (C-VS / in vitro VS)',
       title = 'Cellular vs Inherent Variation Sensitivity') +
  scale_y_continuous(trans = 'log2') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
################################################################################

## Panel 2C/2D: K562 vs HepG2 Fold Change:
################################################################################
# Compute K562/HepG2 fold changes for RBPs with both cell lines:
both_cells$FC_CS_K562_HepG2 = both_cells$IS_K562 / both_cells$IS_HepG2
both_cells$FC_CVS_K562_HepG2 = both_cells$VS_K562 / both_cells$VS_HepG2

## Panel 2C: K562 / HepG2 CS Fold Change:
################################################################################
both_cells_CS = both_cells %>% arrange(desc(FC_CS_K562_HepG2))
both_cells_CS$RBP = factor(both_cells_CS$RBP, levels = both_cells_CS$RBP)

ggplot(both_cells_CS, aes(x = RBP, y = FC_CS_K562_HepG2)) +
  geom_point(fill = 'cornflowerblue', alpha = 1, shape = 21, size = 4) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50') +
  labs(x = 'RBP (ordered by K562/HepG2 CS Fold Change)',
       y = 'Fold Change (K562 CS / HepG2 CS)',
       title = 'K562 vs HepG2 Cellular Specificity') +
  scale_y_continuous(trans = 'log2') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
################################################################################

## Panel 2D: K562 / HepG2 C-VS Fold Change:
################################################################################
both_cells_CVS = both_cells %>% arrange(desc(FC_CVS_K562_HepG2))
both_cells_CVS$RBP = factor(both_cells_CVS$RBP, levels = both_cells_CVS$RBP)

ggplot(both_cells_CVS, aes(x = RBP, y = FC_CVS_K562_HepG2)) +
  geom_point(fill = 'darkgoldenrod1', alpha = 1, shape = 21, size = 4) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50') +
  labs(x = 'RBP (ordered by K562/HepG2 C-VS Fold Change)',
       y = 'Fold Change (K562 C-VS / HepG2 C-VS)',
       title = 'K562 vs HepG2 Cellular Variation Sensitivity') +
  scale_y_continuous(trans = 'log2') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
################################################################################

## Panel 3: CS vs C-VS Correlation for All eCLIP:
################################################################################
# Use representative values (already computed above):
eCLIP_rep = eCLIP_summary %>% filter(!is.na(CS) & !is.na(CVS))

Corr_CS_CVS = cor(eCLIP_rep$CS, eCLIP_rep$CVS, method = 'pearson')

ggplot(eCLIP_rep, aes(x = CS, y = CVS)) +
  geom_point(fill = 'cornflowerblue', alpha = 1, shape = 21, size = 4) +
  geom_smooth(method = 'lm', linetype = 'solid', color = 'cornflowerblue',
              se = FALSE, linewidth = 0.5, alpha = 0.1) +
  labs(x = 'Cellular Specificity',
       y = 'Cellular Variation Sensitivity',
       title = 'CS vs C-VS Correlation (All eCLIP)') +
  # scale_x_continuous() +
  scale_x_continuous(trans = 'log2', limits = c(1, 30)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  annotate('text', x = Inf, y = 0.1, label = sprintf('R = %.2f', Corr_CS_CVS),
           hjust = 1.1, vjust = 0, size = 5, color = 'cornflowerblue')
################################################################################

## Panel 4: Unconventional vs Conventional RBP Boxplots:
################################################################################
# Classify RBPs:
eCLIP_rep$Type = ifelse(eCLIP_rep$RBP %in% unconventional_RBPs,
                        'Unconventional', 'Conventional')

## Panel 4A: CS Boxplot:
################################################################################
ggplot(eCLIP_rep, aes(x = Type, y = CS)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(shape = Type), width = 0.1, size = 3, alpha = 1) +
  labs(x = 'RBP Type',
       y = 'Cellular Specificity',
       title = 'Cellular Specificity: Unconventional vs Conventional') +
  scale_y_continuous(trans = 'log2') +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  geom_signif(comparisons = list(c('Conventional', 'Unconventional')),
              test = 'wilcox.test',
              map_signif_level = FALSE,
              textsize = 5, y_position = max(log2(eCLIP_rep$CS)) + 0.5,
              tip_length = 0.01)
################################################################################

## Panel 4B: C-VS Boxplot:
################################################################################
ggplot(eCLIP_rep, aes(x = Type, y = CVS)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(shape = Type), width = 0.1, size = 3, alpha = 1) +
  labs(x = 'RBP Type',
       y = 'Cellular Variation Sensitivity',
       title = 'Cellular Variation Sensitivity: Unconventional vs Conventional') +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  geom_signif(comparisons = list(c('Conventional', 'Unconventional')),
              test = 'wilcox.test',
              map_signif_level = FALSE,
              textsize = 5, y_position = 0.95, tip_length = 0.01)
################################################################################

## UniProt Domain and Functional Annotation for eCLIP RBPs:
################################################################################
library(httr2)
library(jsonlite)
library(purrr)
library(tidyr)

# Use all eCLIP RBPs:
proteins = unique(eCLIP_summary$RBP)

fetch_uniprot_annotations = function(gene_name) {
  base_url = 'https://rest.uniprot.org/uniprotkb/search'
  query = paste0('(gene:', gene_name, ') AND (organism_id:9606) AND (reviewed:true)')

  resp = request(base_url) %>%
    req_url_query(query = query,
                  fields = paste0('accession,gene_names,protein_name,',
                                  'ft_domain,ft_motif,ft_region,',
                                  'cc_function,cc_subcellular_location,',
                                  'go_p,go_f,go_c,keyword'),
                  format = 'json',
                  size = 1) %>%
    req_perform()

  data = resp_body_json(resp)

  if (length(data$results) == 0) {
    features_df = data.frame(gene = gene_name, accession = NA, type = NA,
                             description = NA, start = NA, end = NA)
    annot_df = data.frame(gene = gene_name, accession = NA, category = NA, value = NA)
    return(list(features = features_df, annotations = annot_df))
  }

  entry = data$results[[1]]
  acc = entry$primaryAccession

  # --- Parse features (domains, motifs, regions) ---
  features = entry$features
  if (length(features) == 0) {
    features_df = data.frame(gene = gene_name, accession = acc, type = NA,
                             description = NA, start = NA, end = NA)
  } else {
    features_df = map_dfr(features, function(f) {
      data.frame(gene = gene_name,
                 accession = acc,
                 type = f$type %||% NA_character_,
                 description = f$description %||% NA_character_,
                 start = f$location$start$value %||% NA_integer_,
                 end = f$location$end$value %||% NA_integer_)
    })
  }

  # --- Parse functional annotations ---
  annot_rows = list()

  # Function (from comments):
  comments = entry$comments
  if (!is.null(comments)) {
    for (cc in comments) {
      if (!is.null(cc$commentType)) {
        if (cc$commentType == 'FUNCTION' && !is.null(cc$texts)) {
          for (txt in cc$texts) {
            annot_rows = c(annot_rows, list(data.frame(
              gene = gene_name, accession = acc,
              category = 'Function', value = txt$value %||% NA_character_)))
          }
        }
        if (cc$commentType == 'SUBCELLULAR LOCATION' && !is.null(cc$subcellularLocations)) {
          for (loc in cc$subcellularLocations) {
            loc_name = loc$location$value %||% NA_character_
            annot_rows = c(annot_rows, list(data.frame(
              gene = gene_name, accession = acc,
              category = 'Subcellular_Location', value = loc_name)))
          }
        }
      }
    }
  }

  # Keywords:
  keywords = entry$keywords
  if (!is.null(keywords)) {
    for (kw in keywords) {
      annot_rows = c(annot_rows, list(data.frame(
        gene = gene_name, accession = acc,
        category = paste0('Keyword_', kw$category %||% 'Unknown'),
        value = kw$name %||% NA_character_)))
    }
  }

  # GO terms:
  go_references = entry$uniProtKBCrossReferences
  if (!is.null(go_references)) {
    for (ref in go_references) {
      if (!is.null(ref$database) && ref$database == 'GO') {
        go_id = ref$id %||% NA_character_
        go_props = ref$properties
        go_term = NA_character_
        if (!is.null(go_props)) {
          for (prop in go_props) {
            if (!is.null(prop$key) && prop$key == 'GoTerm') {
              go_term = prop$value %||% NA_character_
            }
          }
        }
        # Classify by GO aspect:
        go_category = 'GO_Unknown'
        if (!is.na(go_term)) {
          if (startsWith(go_term, 'F:')) { go_category = 'GO_Molecular_Function' }
          else if (startsWith(go_term, 'P:')) { go_category = 'GO_Biological_Process' }
          else if (startsWith(go_term, 'C:')) { go_category = 'GO_Cellular_Component' }
        }
        annot_rows = c(annot_rows, list(data.frame(
          gene = gene_name, accession = acc,
          category = go_category, value = paste0(go_id, ' - ', go_term))))
      }
    }
  }

  if (length(annot_rows) == 0) {
    annot_df = data.frame(gene = gene_name, accession = acc, category = NA, value = NA)
  } else {
    annot_df = do.call(rbind, annot_rows)
  }

  return(list(features = features_df, annotations = annot_df))
}

# Fetch all with delay:
all_features = list()
all_annotations = list()

for (g in proteins) {
  message('Fetching: ', g, ' [', which(proteins == g), '/', length(proteins), ']')
  result = fetch_uniprot_annotations(g)
  all_features[[g]] = result$features
  all_annotations[[g]] = result$annotations
  Sys.sleep(0.5)
}

domain_results = do.call(rbind, all_features)
annotation_results = do.call(rbind, all_annotations)

# Filter to Domain-type annotations:
domains_only = domain_results %>% filter(type == 'Domain')

domains_only %>%
  select(gene, accession, description, start, end) %>%
  print()

# Summarize domain types:
domain_summary = domains_only %>% count(description, sort = TRUE)
print(domain_summary)

# View functional annotation categories:
annotation_results %>% count(category, sort = TRUE) %>% print()

# View GO Molecular Function terms:
annotation_results %>% filter(category == 'GO_Molecular_Function') %>%
  select(gene, value) %>% print()

# View Subcellular Location:
annotation_results %>% filter(category == 'Subcellular_Location') %>%
  select(gene, value) %>% print()

# Keyword Summary 
keyword_summary = annotation_results %>%
  filter(category %in% c('GO_Molecular_Function', 'GO_Biological_Process',
                         'GO_Cellular_Component', 'Keyword_Disease', 'Function')) %>%
  mutate(
    # Clean GO values: "GO:0005737 - C:cytoplasm" -> "cytoplasm"
    value = gsub('^GO:\\d+ - [A-Z]:', '', value),
    # Rename categories:
    category = case_when(
      category == 'GO_Biological_Process' ~ 'Biological_Process',
      category == 'GO_Molecular_Function' ~ 'Molecular_Function',
      category == 'GO_Cellular_Component' ~ 'Cellular_Component',
      category == 'Keyword_Disease' ~ 'Disease',
      category == 'Function' ~ 'Function')) %>%
  group_by(gene, category) %>%
  summarise(value = paste(value, collapse = ', '), .groups = 'drop') %>%
  pivot_wider(id_cols = 'gene', names_from = 'category', values_from = 'value') %>%
  dplyr::rename(RBP = gene) %>%
  dplyr::select(RBP, Cellular_Component, Biological_Process,
                Molecular_Function, Disease, Function)
View(keyword_summary)


################################################################################
