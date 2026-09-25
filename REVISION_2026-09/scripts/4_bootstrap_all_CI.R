################################################################################
## Consolidated CIs for every correlation panel Joe will annotate:
##   Fig 1D  IS vs VS   (n=26 single-RBD RBPs, RBNS)
##   Fig 2C  IS vs CS   (Low / High / All structural-context groups)
##   Fig 2D  VS vs CVS  (Low / High / All)
##   Fig 6A  CS  K562 vs HepG2 (n=76)
##   Fig S6A CVS K562 vs HepG2 (n=76)
## For each: point r; bootstrap 95% percentile CI + % resamples>0 (EXACT enumeration
## for n<=10, else 10,000 Monte-Carlo, seed 1); Fisher-z 95% CI; jackknife most-
## influential RBP (drop-one). Raw scale (matches the quoted figure values).
################################################################################

suppressPackageStartupMessages({ library(readxl); library(readr); library(dplyr); library(RBPSpecificity) })
options(width = 200); set.seed(1)
repo = '/Users/soonyi/Repos/BITS_Specificity/'; dataDir = paste0(repo,'Dataset/Analysis/'); outDir = paste0(repo,'REVISION_2026-09/')
num = function(x) suppressWarnings(as.numeric(x))

## ---- CI machinery ------------------------------------------------------------
wcor = function(x, y, w) { sw = sum(w); mx = sum(w*x)/sw; my = sum(w*y)/sw
  cxy = sum(w*(x-mx)*(y-my)); vx = sum(w*(x-mx)^2); vy = sum(w*(y-my)^2)
  if (vx == 0 || vy == 0) return(NA_real_); cxy/sqrt(vx*vy) }
comps = function(total, parts) {                       # all count vectors summing to total
  if (parts == 1) return(matrix(total,1,1))
  do.call(rbind, lapply(0:total, function(f) cbind(f, comps(total-f, parts-1)))) }
fisher = function(r, n) { z = atanh(r); se = 1/sqrt(n-3); tanh(c(z-1.96*se, z+1.96*se)) }

ci_all = function(x, y, label) {
  ok = is.finite(x) & is.finite(y); x = x[ok]; y = y[ok]; n = length(x); R = cor(x,y)
  if (n <= 10) {                                        # EXACT weighted bootstrap
    C = comps(n, n); logw = lfactorial(n) - rowSums(lfactorial(C)) - n*log(n)
    rr = apply(C, 1, function(c) wcor(x, y, c)); keep = is.finite(rr)
    rr = rr[keep]; w = exp(logw[keep]); w = w/sum(w)
    o = order(rr); rr = rr[o]; w = w[o]; cw = cumsum(w)
    lo = rr[which(cw>=0.025)[1]]; hi = rr[which(cw>=0.975)[1]]; fpos = sum(w[rr>0]); meth = 'exact'
  } else {                                              # Monte-Carlo bootstrap
    bs = replicate(10000, { i = sample.int(n,n,replace=TRUE); suppressWarnings(cor(x[i],y[i])) })
    bs = bs[is.finite(bs)]; ci = quantile(bs, c(.025,.975)); lo=ci[1]; hi=ci[2]; fpos=mean(bs>0); meth='MC-10000'
  }
  fz = fisher(R, n)
  jr = vapply(seq_len(n), function(i) cor(x[-i], y[-i]), numeric(1)); d = jr - R; k = which.max(abs(d))
  data.frame(panel = label, n = n, r = round(R,3),
             boot_lo = round(lo,3), boot_hi = round(hi,3), frac_pos = round(fpos,4), boot_method = meth,
             fisher_lo = round(fz[1],3), fisher_hi = round(fz[2],3),
             jack_worst = names(x)[k] %||% k, jack_r_without = round(jr[k],3), jack_delta = round(d[k],3))
}
`%||%` = function(a,b) if (is.null(a)) b else a

## ---- Fig 1D : IS vs VS, 26 single-RBD RBPs -----------------------------------
rb = read_csv(paste0(dataDir,'RBNS/RBNS_singleRBD_normalized_5mer.csv'), show_col_types=FALSE)
rbps = colnames(rb)[-1]
isvs = sapply(rbps, function(p){ t = data.frame(MOTIF=rb$Motif, Score=rb[[p]])
  top = t$MOTIF[which.max(t$Score)]                        # VS = returnMS on RBNS (0.99.0 has no returnVS)
  c(IS = suppressMessages(returnIS(t, top)), VS = suppressMessages(returnMS(t, top, output_type='number'))) })
IS = setNames(isvs['IS',], rbps); VS = setNames(isvs['VS',], rbps)
fig1d = ci_all(IS, VS, 'Fig1D IS-VS (n=26)')

## ---- Fig 2C / 2D : from Table S5 + structural context ------------------------
s5 = read_excel(paste0(repo,'Dataset/Table_S5.xlsx'))
s5 = data.frame(RBP=s5$RBP, IS=num(s5$IS), VS=num(s5$VS), CS=num(s5$CS), CVS=num(s5$CVS))
st = read.csv(paste0(dataDir,'sample_table_eCLIP.csv'), fileEncoding='UTF-8-BOM'); st$RBP0 = sub('_[0-9]+$','',st$RBP)
struct = st %>% distinct(RBP0, structure_context) %>% rename(RBP=RBP0) %>%
  filter(structure_context %in% c('H','L')) %>% mutate(group=ifelse(structure_context=='H','High','Low'))
f2 = struct %>% inner_join(s5, by='RBP')
named = function(v, nm){ names(v) = nm; v }
g = function(grp) if (grp=='All') f2 else f2 %>% filter(group==grp)
mk = function(panel, xcol, ycol) bind_rows(lapply(c('Low','High','All'), function(grp){
  d = g(grp); ci_all(named(d[[xcol]], d$RBP), named(d[[ycol]], d$RBP), sprintf('%s %s', panel, grp)) }))
fig2c = mk('Fig2C IS-CS', 'IS', 'CS')
fig2d = mk('Fig2D VS-CVS', 'VS', 'CVS')

## ---- Fig 6A / S6A : K562 vs HepG2 (n=76) from Task A -------------------------
a = read_csv(paste0(outDir,'4_bootstrap_cellline_reproducibility.csv'), show_col_types=FALSE) %>%
  filter(record_type=='per_RBP', set=='Fig6_set')
fig6a  = ci_all(named(a$CS_K562, a$RBP),  named(a$CS_HepG2, a$RBP),  'Fig6A CS K562-HepG2 (n=76)')
figs6a = ci_all(named(a$CVS_K562, a$RBP), named(a$CVS_HepG2, a$RBP), 'FigS6A CVS K562-HepG2 (n=76)')

out = bind_rows(fig1d, fig2c, fig2d, fig6a, figs6a)
write_csv(out, paste0(outDir,'4_bootstrap_all_CI.csv'), na='')
print(out, row.names=FALSE)
