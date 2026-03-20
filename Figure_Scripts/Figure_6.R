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
library(ggrepel)

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

# Merge RBNS IS for coloring:
both_cells = merge(both_cells, RBNS_metrics[, c('RBP', 'IS')], by = 'RBP', all.x = TRUE)
colnames(both_cells)[colnames(both_cells) == 'IS'] = 'RBNS_IS'

trend_model = lm(IS_K562 ~ IS_HepG2, data = both_cells)
both_cells$Trend_Residual = rstandard(trend_model)
both_cells$Trend_Outlier = abs(both_cells$Trend_Residual) > 2

Corr_CS = cor(both_cells$IS_K562, both_cells$IS_HepG2, method = 'pearson')

ggplot(both_cells, aes(x = IS_K562, y = IS_HepG2, fill = (RBNS_IS))) +
  geom_point(alpha = 1, shape = 21, size = 4) +
  scale_fill_viridis_c(name = 'RBNS IS', na.value = 'grey70', option = 'viridis') +
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
# Merge RBNS VS for coloring:
both_cells = merge(both_cells, RBNS_metrics[, c('RBP', 'VS')], by = 'RBP', all.x = TRUE)
colnames(both_cells)[colnames(both_cells) == 'VS'] = 'RBNS_VS'

trend_model = lm(VS_K562 ~ VS_HepG2, data = both_cells)
both_cells$Trend_Residual = rstandard(trend_model)
both_cells$Trend_Outlier = abs(both_cells$Trend_Residual) > 2

Corr_CVS = cor(both_cells$VS_K562, both_cells$VS_HepG2, method = 'pearson')

ggplot(both_cells, aes(x = VS_K562, y = VS_HepG2, fill = RBNS_VS)) +
  geom_point(alpha = 1, shape = 21, size = 4) +
  scale_fill_viridis_c(name = 'RBNS VS', na.value = 'grey70', option = 'viridis') +
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

# Identify trendline outliers using standard residuals
trend_model <- lm(CVS ~ log2(CS), data = eCLIP_rep)
eCLIP_rep$Trend_Residual <- rstandard(trend_model)
eCLIP_rep$Trend_Outlier <- abs(eCLIP_rep$Trend_Residual) > 2

Corr_CS_CVS = cor(eCLIP_rep$CS, eCLIP_rep$CVS, method = 'pearson')

ggplot(eCLIP_rep, aes(x = CS, y = CVS)) +
  geom_point(aes(fill = Trend_Outlier), alpha = 1, shape = 21, size = 4) +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  geom_smooth(method = 'lm', linetype = 'solid', color = 'grey',
              se = FALSE, linewidth = 0.5, alpha = 0.1) +
  geom_text_repel(data = subset(eCLIP_rep, Trend_Outlier), aes(label = RBP), 
                  size = 3, box.padding = 0.8, min.segment.length = 0, segment.color = 'grey50', max.overlaps = 50) +
  labs(x = 'Cellular Specificity',
       y = 'Cellular Variation Sensitivity',
       title = 'CS vs C-VS Correlation (All eCLIP)') +
  scale_x_continuous(trans = 'log2', limits = c(1, 30)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  guides(fill = "none") +
  annotate('text', x = Inf, y = 0.1, label = sprintf('R = %.2f', Corr_CS_CVS),
           hjust = 1.1, vjust = 0, size = 5, color = 'black')
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
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  geom_signif(comparisons = list(c('Conventional', 'Unconventional')),
              test = 'wilcox.test',
              map_signif_level = FALSE,
              textsize = 5, y_position = 0.95, tip_length = 0.01)
################################################################################

## Panel 5: SR vs HNRNP IS and VS Comparison:
################################################################################
SR_proteins = c('SRSF2', 'SRSF4', 'SRSF5', 'SRSF8', 'SRSF9', 'SRSF10', 'SRSF11')
HNRNP_proteins = c('HNRNPA0', 'HNRNPA2B1', 'HNRNPC', 'HNRNPCL1', 'HNRNPD0',
                   'HNRNPDL', 'HNRNPF', 'HNRNPH2', 'HNRNPK', 'HNRNPL')

SR_HNRNP = RBNS_metrics %>%
  filter(RBP %in% c(SR_proteins, HNRNP_proteins)) %>%
  mutate(Family = ifelse(RBP %in% SR_proteins, 'SR', 'HNRNP'))

## IS Boxplot:
pos_jitter = position_jitter(width = 0.1, seed = 42)
ggplot(SR_HNRNP, aes(x = Family, y = IS, label = RBP)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter, size = 3, alpha = 1) +
  geom_text_repel(position = pos_jitter, size = 3, max.overlaps = 20) +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  labs(x = 'Protein Family',
       y = 'Inherent Specificity',
       title = 'IS: SR vs HNRNP (RBNS 5-mer)') +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  geom_signif(comparisons = list(c('HNRNP', 'SR')),
              test = 'wilcox.test',
              map_signif_level = FALSE,
              textsize = 5, tip_length = 0.01)

## VS Boxplot:
pos_jitter = position_jitter(width = 0.1, seed = 42)
ggplot(SR_HNRNP, aes(x = Family, y = VS, label = RBP)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter, size = 3, alpha = 1) +
  geom_text_repel(position = pos_jitter, size = 3, max.overlaps = 20) +
  labs(x = 'Protein Family',
       y = 'Variation Sensitivity',
       title = 'VS: SR vs HNRNP (RBNS 5-mer)') +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14)) +
  geom_signif(comparisons = list(c('HNRNP', 'SR')),
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

# Summarize domain types:
domain_summary = domains_only %>% count(description, sort = TRUE)

# Per-RBP domain table:
domain_per_rbp = domains_only %>%
  group_by(gene) %>%
  summarise(Domains = paste(description, collapse = ', '), .groups = 'drop') %>%
  dplyr::rename(RBP = gene)

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
################################################################################

## Meta Summary Table:
################################################################################
# Load RBNS enrichment data for top motifs:
RBNS_enrich = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'),
                       col_names = TRUE, show_col_types = FALSE)
RBNS_enrich = data.frame(RBNS_enrich)
RBNS_RBPs = colnames(RBNS_enrich)[2:ncol(RBNS_enrich)]

# Helper: get top N motifs as comma-separated string:
get_top_motifs = function(motifs, scores, N = 10) {
  idx = order(scores, decreasing = TRUE)[1:min(N, length(scores))]
  paste(motifs[idx], collapse = ', ')
}

# Get RBNS top 10 motifs per RBP:
rbns_top10 = data.frame(RBP = RBNS_RBPs, Top10_Motif_RBNS = NA_character_,
                        stringsAsFactors = FALSE)
for (i in seq_along(RBNS_RBPs)) {
  rbp = RBNS_RBPs[i]
  rbns_top10$Top10_Motif_RBNS[i] = get_top_motifs(RBNS_enrich$Motif,
                                                    RBNS_enrich[[rbp]])
}

# Get eCLIP top 10 motifs per RBP (use averaged enrichment if both cell lines):
eCLIP_RBPs = unique(eCLIP_summary$RBP)
eCLIP_top10 = data.frame(RBP = eCLIP_RBPs, Top10_Motif_eCLIP = NA_character_,
                         stringsAsFactors = FALSE)
eCLIP_enrichDir = paste0(baseDir, 'eCLIP_all/output/', DATE, '/')
extension = 25

for (i in seq_along(eCLIP_RBPs)) {
  rbp = eCLIP_RBPs[i]
  has_K562 = !is.na(eCLIP_summary$IS_K562[eCLIP_summary$RBP == rbp])
  has_HepG2 = !is.na(eCLIP_summary$IS_HepG2[eCLIP_summary$RBP == rbp])

  enrich_data = NULL

  if (has_K562 && has_HepG2) {
    # Average both cell lines:
    K562_file = paste0(eCLIP_enrichDir, 'K562_eCLIP_Enrichment_', rbp, '_', K, 'mer.csv')
    HepG2_file = paste0(eCLIP_enrichDir, 'HepG2_eCLIP_Enrichment_', rbp, '_', K, 'mer.csv')
    if (file.exists(K562_file) && file.exists(HepG2_file)) {
      K562_data = read_csv(K562_file, col_names = TRUE, show_col_types = FALSE)
      HepG2_data = read_csv(HepG2_file, col_names = TRUE, show_col_types = FALSE)
      K562_data = K562_data %>% arrange(MOTIF)
      HepG2_data = HepG2_data %>% arrange(MOTIF)
      enrich_data = data.frame(MOTIF = K562_data$MOTIF,
                               Score = (K562_data$Score + HepG2_data$Score) / 2)
    }
  } else if (has_K562) {
    K562_file = paste0(eCLIP_enrichDir, 'K562_eCLIP_Enrichment_', rbp, '_', K, 'mer.csv')
    if (file.exists(K562_file)) {
      enrich_data = read_csv(K562_file, col_names = TRUE, show_col_types = FALSE)
    }
  } else if (has_HepG2) {
    HepG2_file = paste0(eCLIP_enrichDir, 'HepG2_eCLIP_Enrichment_', rbp, '_', K, 'mer.csv')
    if (file.exists(HepG2_file)) {
      enrich_data = read_csv(HepG2_file, col_names = TRUE, show_col_types = FALSE)
    }
  }

  if (!is.null(enrich_data)) {
    enrich_data$MOTIF = gsub('T', 'U', enrich_data$MOTIF)
    eCLIP_top10$Top10_Motif_eCLIP[i] = get_top_motifs(enrich_data$MOTIF,
                                                       enrich_data$Score)
  }
}

# Compile meta table:
# Load RBP Domain Summary:
rbd_info = read_csv(paste0(baseDir, 'eCLIP_all/output/eCLIP_RBPs_Uniprot_RBD_Table.csv'),
                    col_names = TRUE, show_col_types = FALSE)
rbd_info = data.frame(rbd_info)
rbd_info[] = lapply(rbd_info, function(x) {
  if (is.character(x)) gsub('\u00A0', ' ', x) else x
})
# Rename columns for easy use in R:
colnames(rbd_info) = c('RBP', 'Domains', 'Domain_Architecture', 'Primary_RBD', 'Other_Domains', 'Disordered_Regions')

meta_table = eCLIP_summary %>%
  dplyr::select(RBP, CS, CVS) %>%
  merge(RBNS_metrics, by = 'RBP', all = TRUE) %>%
  merge(rbns_top10, by = 'RBP', all.x = TRUE) %>%
  merge(eCLIP_top10, by = 'RBP', all.x = TRUE) %>%
  merge(rbd_info, by = 'RBP', all.x = TRUE) %>%
  mutate(idx = match(RBP, eCLIP_summary$RBP),
         RBNS = ifelse(!is.na(IS), 'Yes', 'No'),
         eCLIP = ifelse(!is.na(CS), 'Yes', 'No'),
         K562 = ifelse(!is.na(idx) & !is.na(eCLIP_summary$IS_K562[idx]),
                       'Yes', 'No'),
         HepG2 = ifelse(!is.na(idx) & !is.na(eCLIP_summary$IS_HepG2[idx]),
                        'Yes', 'No')) %>%
  dplyr::select(RBP, RBNS, eCLIP, IS, VS, Top10_Motif_RBNS, CS, CVS, Top10_Motif_eCLIP, K562, HepG2,
                Domains, Domain_Architecture, Primary_RBD, Other_Domains, Disordered_Regions)

write.csv(meta_table, paste0(baseDir, '5mer_MetaAnalysis_Output.csv'), row.names = FALSE, na = '')
################################################################################

## Panel 6: CS and CVS Boxplots by RBP Domain Categorizations
################################################################################
# Merge domain info into eCLIP_rep if not already present
if (!'Primary_RBD' %in% colnames(eCLIP_rep)) {
  eCLIP_rep = merge(eCLIP_rep, rbd_info, by = 'RBP', all.x = TRUE)
}

# 1) Primary classification
eCLIP_rep = eCLIP_rep %>%
  mutate(Primary_Classification = case_when(
    RBP %in% c('IGF2BP1', 'IGF2BP2', 'IGF2BP3') ~ 'RRM',
    RBP == 'DDX43' ~ 'Helicase',
    grepl('RRM', Primary_RBD) ~ 'RRM',
    grepl('KH', Primary_RBD) ~ 'KH',
    grepl('Helicase', Primary_RBD) ~ 'Helicase',
    grepl('dsRBM', Primary_RBD) ~ 'dsRBD',
    TRUE ~ ifelse(!is.na(Disordered_Regions) & grepl('YES', Disordered_Regions), 
                  'Others (Disordered)', 'Others (No Disordered Regions)')
  ))

# Ensure factor order for plotting x-axis
eCLIP_rep$Primary_Classification = factor(eCLIP_rep$Primary_Classification, 
                                          levels = c('RRM', 'KH', 'Helicase', 'dsRBD', 
                                                     'Others (Disordered)', 'Others (No Disordered Regions)'))

# 2) Secondary classification (Shape)
eCLIP_rep = eCLIP_rep %>%
  mutate(Shape_Class = case_when(
    grepl('Others', Primary_Classification) ~ 'Others',
    !is.na(Domain_Architecture) & Domain_Architecture == 'Single' ~ 'Single',
    !is.na(Domain_Architecture) & Domain_Architecture == 'Multi' ~ 'Multi',
    !is.na(Domain_Architecture) & Domain_Architecture == 'Mixed' ~ 'Mixed',
    TRUE ~ 'Others'
  ))

eCLIP_rep$Shape_Class = factor(eCLIP_rep$Shape_Class, levels = c('Single', 'Multi', 'Mixed', 'Others'))

# Define shapes: 16 = circle, 15 = square(rectangle), 17 = triangle, 18 = diamond
shape_mapping = c('Single' = 16, 'Multi' = 15, 'Mixed' = 17, 'Others' = 18)

## Panel 6A: CS Boxplot by Domain
kw_cs = kruskal.test(CS ~ Primary_Classification, data = eCLIP_rep)
pos_dodge_jitter = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75, seed = 42)
ggplot(eCLIP_rep, aes(x = Primary_Classification, y = CS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(data = subset(eCLIP_rep, !grepl('Others', Primary_Classification)),
             aes(shape = Shape_Class, fill = Shape_Class, group = Shape_Class), 
             position = pos_dodge_jitter, size = 3, alpha = 0.7) +
  geom_point(data = subset(eCLIP_rep, grepl('Others', Primary_Classification)),
             aes(shape = Shape_Class, fill = Shape_Class), 
             position = position_jitter(width = 0.2, seed = 42), size = 3, alpha = 0.7) +
  scale_shape_manual(values = shape_mapping) +
  labs(x = 'Primary RBD Classification',
       y = 'Cellular Specificity',
       title = 'CS: Domain Categorizations',
       shape = 'Domain Architecture') +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12)) +
  annotate('text', x = Inf, y = 64, label = sprintf('K-W p = %.2e', kw_cs$p.value),
           hjust = 1.1, vjust = 1.5, size = 5, fontface = 'bold')

## Panel 6B: CVS Boxplot by Domain
kw_cvs = kruskal.test(CVS ~ Primary_Classification, data = eCLIP_rep)
ggplot(eCLIP_rep, aes(x = Primary_Classification, y = CVS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(data = subset(eCLIP_rep, !grepl('Others', Primary_Classification)),
             aes(shape = Shape_Class, fill = Shape_Class, group = Shape_Class), 
             position = pos_dodge_jitter, size = 3, alpha = 0.7) +
  geom_point(data = subset(eCLIP_rep, grepl('Others', Primary_Classification)),
             aes(shape = Shape_Class, fill = Shape_Class), 
             position = position_jitter(width = 0.2, seed = 42), size = 3, alpha = 0.7) +
  scale_shape_manual(values = shape_mapping) +
  labs(x = 'Primary RBD Classification',
       y = 'Cellular Variation Sensitivity',
       title = 'C-VS: Domain Categorizations',
       shape = 'Domain Architecture') +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 12)) +
  annotate('text', x = Inf, y = 1.0, label = sprintf('K-W p = %.2e', kw_cvs$p.value),
           hjust = 1.1, vjust = 1.5, size = 5, fontface = 'bold')
################################################################################

## Panel 7: CS and CVS Boxplots by Disordered Regions (All RBPs)
################################################################################
eCLIP_rep = eCLIP_rep %>%
  mutate(Disordered_All = ifelse(!is.na(Disordered_Regions) & grepl('YES', Disordered_Regions), 
                                 'Disordered', 'No Disordered Regions'))

eCLIP_rep$Disordered_All = factor(eCLIP_rep$Disordered_All, levels = c('No Disordered Regions', 'Disordered'))

## Panel 7A: CS Boxplot by Disordered Regions
ggplot(eCLIP_rep, aes(x = Disordered_All, y = CS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  labs(x = 'Disordered Regions',
       y = 'Cellular Specificity',
       title = 'CS: Disordered vs No Disordered Regions (All RBPs)') +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold')) +
  geom_signif(comparisons = list(c('No Disordered Regions', 'Disordered')),
              test = 'wilcox.test', map_signif_level = FALSE,
              textsize = 5, tip_length = 0.01)

## Panel 7B: CVS Boxplot by Disordered Regions
ggplot(eCLIP_rep, aes(x = Disordered_All, y = CVS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  labs(x = 'Disordered Regions',
       y = 'Cellular Variation Sensitivity',
       title = 'C-VS: Disordered vs No Disordered Regions (All RBPs)') +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold')) +
  geom_signif(comparisons = list(c('No Disordered Regions', 'Disordered')),
              test = 'wilcox.test', map_signif_level = FALSE,
              textsize = 5, y_position = 0.95, tip_length = 0.01)
################################################################################

## New Panels (Graph C, D, E): Domain Comparisons with Outlier Labeling
################################################################################
# Helper function to find boxplot outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm=TRUE) - 1.5 * IQR(x, na.rm=TRUE) | 
         x > quantile(x, 0.75, na.rm=TRUE) + 1.5 * IQR(x, na.rm=TRUE))
}

eCLIP_rep = eCLIP_rep %>%
  mutate(
    # Graph C: Primary RBD (5 Categories)
    Primary_Class = case_when(
      RBP %in% c('IGF2BP1', 'IGF2BP2', 'IGF2BP3') ~ 'RRM',
      RBP == 'DDX43' ~ 'Helicase',
      grepl('RRM', Primary_RBD) ~ 'RRM',
      grepl('KH', Primary_RBD) ~ 'KH',
      grepl('Helicase', Primary_RBD) ~ 'Helicase',
      grepl('dsRBM', Primary_RBD) ~ 'dsRBD',
      TRUE ~ 'Others'
    ),
    # Graph D: Domain Architecture (3 Categories)
    Arch_Class = case_when(
      Domain_Architecture %in% c('Single', 'Multi', 'Mixed') ~ Domain_Architecture,
      TRUE ~ NA_character_
    ),
    # Graph E: Disordered Regions (2 Categories)
    Disordered_Class = ifelse(!is.na(Disordered_Regions) & grepl('YES', Disordered_Regions), 
                              'Disordered', 'No Disordered Regions')
  ) %>%
  # Calculate outliers per group for labeling
  group_by(Primary_Class) %>%
  mutate(Outlier_C_CS = is_outlier(log2(CS)), Outlier_C_CVS = is_outlier(CVS)) %>%
  group_by(Arch_Class) %>%
  mutate(Outlier_D_CS = is_outlier(log2(CS)), Outlier_D_CVS = is_outlier(CVS)) %>%
  group_by(Disordered_Class) %>%
  mutate(Outlier_E_CS = is_outlier(log2(CS)), Outlier_E_CVS = is_outlier(CVS)) %>%
  ungroup()

eCLIP_rep$Primary_Class = factor(eCLIP_rep$Primary_Class, levels = c('RRM', 'KH', 'Helicase', 'dsRBD', 'Others'))
eCLIP_rep$Arch_Class = factor(eCLIP_rep$Arch_Class, levels = c('Single', 'Multi', 'Mixed'))
eCLIP_rep$Disordered_Class = factor(eCLIP_rep$Disordered_Class, levels = c('No Disordered Regions', 'Disordered'))

## Graph C: Primary RBD (5 bars)
pos_jitter_new = position_jitter(width = 0.2, seed = 42)

ggplot(eCLIP_rep, aes(x = Primary_Class, y = CS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter_new, size = 3, alpha = 0.7) +
  geom_text_repel(data = eCLIP_rep, aes(label = ifelse(Outlier_C_CS, RBP, NA_character_)), position = pos_jitter_new, size = 3, max.overlaps = 50,
                  min.segment.length = 0, box.padding = 0.8, segment.color = 'grey50') +
  labs(x = 'Primary RBD', y = 'Cellular Specificity', title = 'Graph C: Primary RBD (CS)') +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'))

ggplot(eCLIP_rep, aes(x = Primary_Class, y = CVS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter_new, size = 3, alpha = 0.7) +
  geom_text_repel(data = eCLIP_rep, aes(label = ifelse(Outlier_C_CVS, RBP, NA_character_)), position = pos_jitter_new, size = 3, max.overlaps = 50,
                  min.segment.length = 0, box.padding = 0.8, segment.color = 'grey50') +
  labs(x = 'Primary RBD', y = 'Cellular Variation Sensitivity', title = 'Graph C: Primary RBD (C-VS)') +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'))

## Graph D: Domain Architecture (3 bars)
d_data = subset(eCLIP_rep, !is.na(Arch_Class))

ggplot(d_data, aes(x = Arch_Class, y = CS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter_new, size = 3, alpha = 0.7) +
  geom_text_repel(data = d_data, aes(label = ifelse(Outlier_D_CS, RBP, NA_character_)), position = pos_jitter_new, size = 3, max.overlaps = 50,
                  min.segment.length = 0, box.padding = 0.8, segment.color = 'grey50') +
  labs(x = 'Domain Architecture', y = 'Cellular Specificity', title = 'Graph D: Domain Architecture (CS)') +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'))

ggplot(d_data, aes(x = Arch_Class, y = CVS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter_new, size = 3, alpha = 0.7) +
  geom_text_repel(data = d_data, aes(label = ifelse(Outlier_D_CVS, RBP, NA_character_)), position = pos_jitter_new, size = 3, max.overlaps = 50,
                  min.segment.length = 0, box.padding = 0.8, segment.color = 'grey50') +
  labs(x = 'Domain Architecture', y = 'Cellular Variation Sensitivity', title = 'Graph D: Domain Architecture (C-VS)') +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'))

## Graph E: Disordered Regions (2 bars)
ggplot(eCLIP_rep, aes(x = Disordered_Class, y = CS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter_new, size = 3, alpha = 0.7) +
  geom_text_repel(data = eCLIP_rep, aes(label = ifelse(Outlier_E_CS, RBP, NA_character_)), position = pos_jitter_new, size = 3, max.overlaps = 50,
                  min.segment.length = 0, box.padding = 0.8, segment.color = 'grey50') +
  labs(x = 'Disordered Regions', y = 'Cellular Specificity', title = 'Graph E: Disordered Regions (CS)') +
  scale_y_continuous(trans = 'log2', limits = c(1, 64)) +
  theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'))

ggplot(eCLIP_rep, aes(x = Disordered_Class, y = CVS)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_point(position = pos_jitter_new, size = 3, alpha = 0.7) +
  geom_text_repel(data = eCLIP_rep, aes(label = ifelse(Outlier_E_CVS, RBP, NA_character_)), position = pos_jitter_new, size = 3, max.overlaps = 50,
                  min.segment.length = 0, box.padding = 0.8, segment.color = 'grey50') +
  labs(x = 'Disordered Regions', y = 'Cellular Variation Sensitivity', title = 'Graph E: Disordered Regions (C-VS)') +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'))
################################################################################
