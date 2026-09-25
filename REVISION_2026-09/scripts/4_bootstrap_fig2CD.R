################################################################################
## Task G: bootstrap and jackknife of the Fig. 2C and Fig. 2D correlations
## (Reviewer 2, comment 5)
##
## Panels are split into High/Low structural-context groups (9 + 9 RBPs). The manuscript's
## quoted positive correlations are the LOW group (Fig 2C IS-CS r=0.94; Fig 2D VS-CVS r=0.84),
## so that is the primary target; High and the pooled 18-RBP set are reported alongside.
## Fig 2C uses RAW IS vs CS (matches the reported r); Fig 2D uses VS vs CVS (linear).
## 10,000 bootstrap resamples with replacement over RBPs. seed = 1.
################################################################################

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(readr) })
options(width = 200)
set.seed(1)

repo    = '/Users/soonyi/Repos/BITS_Specificity/'
dataDir = paste0(repo, 'Dataset/Analysis/')
outDir  = paste0(repo, 'REVISION_2026-09/')
num = function(x) suppressWarnings(as.numeric(x))
NBOOT = 10000

s5 = read_excel(paste0(repo, 'Dataset/Table_S5.xlsx'))
s5 = data.frame(RBP = s5$RBP, IS = num(s5$IS), VS = num(s5$VS), CS = num(s5$CS), CVS = num(s5$CVS))
st = read.csv(paste0(dataDir, 'sample_table_eCLIP.csv'), fileEncoding = 'UTF-8-BOM')
st$RBP0 = sub('_[0-9]+$', '', st$RBP)
struct = st %>% distinct(RBP0, structure_context) %>% rename(RBP = RBP0) %>%
  filter(structure_context %in% c('H', 'L')) %>% mutate(group = ifelse(structure_context == 'H', 'High', 'Low'))
fig2 = struct %>% inner_join(s5, by = 'RBP')

## Bootstrap + jackknife for one (panel, group, x, y) set
################################################################################
boot_jack = function(panel, grp, xcol, ycol, df) {
  d = if (grp == 'All') df else df %>% filter(group == grp)
  x = d[[xcol]]; y = d[[ycol]]; n = length(x)
  R = cor(x, y)
  # bootstrap
  bs = numeric(NBOOT)
  for (i in seq_len(NBOOT)) {
    idx = sample.int(n, n, replace = TRUE)
    bs[i] = suppressWarnings(cor(x[idx], y[idx]))
  }
  bs = bs[is.finite(bs)]
  ci = quantile(bs, c(0.025, 0.975))
  # jackknife
  jr = vapply(seq_len(n), function(i) cor(x[-i], y[-i]), numeric(1))
  delta = jr - R
  ord = order(abs(delta), decreasing = TRUE)[1:3]
  summ = data.frame(record_type = 'bootstrap', panel = panel, group = grp, x = xcol, y = ycol, n = n,
                    point_R = round(R, 4), boot_mean = round(mean(bs), 4), boot_SE = round(sd(bs), 4),
                    ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4),
                    n_valid_boot = length(bs), frac_boot_gt0 = round(mean(bs > 0), 4))
  jackrows = data.frame(record_type = 'jackknife', panel = panel, group = grp, x = xcol, y = ycol, n = n,
                        point_R = round(R, 4),
                        dropped_RBP = d$RBP[ord], R_without = round(jr[ord], 4),
                        deltaR = round(delta[ord], 4), rank = 1:3)
  list(summ = summ, jack = jackrows)
}

sets = list(
  list('Fig2C', 'Low',  'IS',  'CS'),   # primary: quoted r = 0.94
  list('Fig2C', 'High', 'IS',  'CS'),
  list('Fig2C', 'All',  'IS',  'CS'),
  list('Fig2D', 'Low',  'VS',  'CVS'),  # primary: quoted r = 0.84
  list('Fig2D', 'High', 'VS',  'CVS'),
  list('Fig2D', 'All',  'VS',  'CVS')
)
boot_all = list(); jack_all = list()
for (s in sets) {
  r = boot_jack(s[[1]], s[[2]], s[[3]], s[[4]], fig2)
  boot_all[[length(boot_all)+1]] = r$summ
  jack_all[[length(jack_all)+1]] = r$jack
}
boot_tab = bind_rows(boot_all); jack_tab = bind_rows(jack_all)

res = bind_rows(boot_tab, jack_tab)
write_csv(res, paste0(outDir, '4_bootstrap_fig2CD.csv'), na = '')

cat('=== Bootstrap (', NBOOT, 'resamples, seed 1). Fig2C = raw IS vs CS; Fig2D = VS vs CVS ===\n')
print(boot_tab)
cat('\n=== Jackknife: 3 RBPs whose removal changes R most (per panel/group) ===\n')
print(jack_tab)
