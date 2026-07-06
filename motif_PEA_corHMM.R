library(corHMM)
library(ape)    # for read.tree

# 1. Read in tree and trait as you already have:
tree       <- read.tree("/Users/user/Documents/PlantProject/Motif_gene_log2FC_data_v1/BT_tree_id_pruned_renamed.tre")
trait_data <- read.csv(
  "/Users/user/Documents/PlantProject/Motif_gene_log2FC_data_v1/trait_values_old_pruned_renamed.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)
thr = 2

# 2. Binarize HSF at threshold = 1
binary_states <- ifelse(trait_data$HSFA > thr, 2, 1)

# 3. Prepare a two‐column data.frame for corHMM:
#    First column must match your tip labels; second is the integer state.
cor_data <- data.frame(
  Species = rownames(trait_data),
  State   = binary_states,
  row.names = NULL
)

# 4. Run ancestral‐state reconstruction.

rc = 1
#    Here’s an example under the “ARD” (all‐rates‐different) model:
res_ARD <- corHMM(
  phy   = tree,
  data  = cor_data,
  rate.cat = rc,
  model = "ARD"
)

# ER model with 1 rate‐category
res_ER  <- corHMM(
  phy      = tree,
  data     = cor_data,
  rate.cat = rc,
  model    = "ER"
)


# 5. Inspect / plot results
print(res_ARD$rates)        # estimated transition rates
print(res_ARD$lik.anc)      # likelihoods at nodes

# Plot with pie‐charts at each internal node:
plotRECON(
  phy         = res_ER$phy,
  likelihoods = res_ER$states,
  pie.cex     = 0.6,
  show.tip.label = TRUE
)

library(ape)  # for plot.phylo() & nodelabels()

# ---- Plot ARD
# 1) define a 2‐color scheme (state 1 = “low”, state 2 = “high”)
piecols  <- c("lightblue","salmon")
pie_cex  <- 0.6

# 2) plot the bare tree
plot(
  res_ARD$phy,
  show.tip.label = TRUE,
  cex            = 0.7
)

# 3) add marginal‐probability pies at each internal node
nodelabels(
  pie    = res_ARD$states,  # your (Nnode × 2) matrix
  piecol = piecols,
  cex    = pie_cex,
  frame  = "none"
)

# 4) add a bottom‐left legend that won’t block the tree
legend(
  "bottomleft",
  legend = c(paste0("HSFA ≤ ", thr),paste0("HSFA > ", thr)),
  fill   = piecols,
  border = NA,
  bty    = "n",
  cex    = 0.8
)
 # ---- Plot ER
# 1) define a 2‐color scheme (state 1 = “low”, state 2 = “high”)
piecols  <- c("lightblue","salmon")
pie_cex  <- 0.6

# 2) plot the bare tree
plot(
  res_ER$phy,
  show.tip.label = TRUE,
  cex            = 0.7
)

# 3) add marginal‐probability pies at each internal node
nodelabels(
  pie    = res_ER$states,  # your (Nnode × 2) matrix
  piecol = piecols,
  cex    = pie_cex,
  frame  = "none"
)

# 4) add a bottom‐left legend that won’t block the tree
legend(
  "bottomleft",
  legend = c(paste0("HSFA ≤ ", thr),paste0("HSFA > ", thr)),
  fill   = piecols,
  border = NA,
  bty    = "n",
  cex    = 0.8
)

# ---- Displace rates

# ER model has a single rate:
rates_ER <- res_ER$solution
cat("ER model fitted rate (q₁₂ = q₂₁):", rates_ER, "\n")

# ARD model has two rates:
Qmat <- res_ARD$solution
q12  <- Qmat[1,2]   # state 1 → 2
q21  <- Qmat[2,1]   # state 2 → 1

cat(" q₁₂ =", q12, "\n",
    " q₂₁ =", q21, "\n")



# ---- Compare AIC
# 1. Extract AICc values
aicc_ER  <- res_ER$AICc
aicc_ARD <- res_ARD$AICc

# 2. Compute ΔAICc = AICc(ER) − AICc(ARD)
deltaAICc <- aicc_ER - aicc_ARD

cat("AICc ER model: ",  aicc_ER,  "\n",
    "AICc ARD model:", aicc_ARD, "\n",
    "ΔAICc (ER − ARD):", deltaAICc, "\n")
# 3. Compute Akaike weights
aiccs    <- c(ER = aicc_ER, ARD = aicc_ARD)
w        <- exp(-0.5 * (aiccs - min(aiccs)))
weights  <- w / sum(w)

print(weights)

# ----Plot info on tree
library(ape)

# 1. Extract your info vector and tree
info_vals <- res_ARD$states.info
phy       <- res_ARD$phy

# 2. Build a 100‐step blue→white→red palette
pal <- colorRampPalette(c("blue","white","red"))(100)

# 3. Map each info value into one of the 100 colors
#    (we cut the range of info_vals into 100 bins)
col_idx     <- as.numeric(cut(info_vals,
                              breaks = 100,
                              include.lowest = TRUE))
node_colors <- pal[col_idx]

# 4. Plot the bare tree
plot(phy, show.tip.label=TRUE, cex=0.7)

# 5. Add a colored circle at each internal node
#    ape::nodelabels by default labels nodes Ntip+1 … Ntip+Nnode
nodelabels(
  pch    = 21,
  bg     = node_colors,
  cex    = 1.2,
  frame  = "none"
)

# 6. Add a continuous color bar legend in the bottom‐left
#    We'll draw a little horizontal strip of colored rectangles
#    and annotate it with five tick‐values.

# figure out user coords
usr <- par("usr")
x_left  <- usr[1] + 0.05*(usr[2]-usr[1])
x_right <- x_left + 0.2*(usr[2]-usr[1])
y_bot   <- usr[3] + 0.05*(usr[4]-usr[3])
y_top   <- y_bot   + 0.05*(usr[4]-usr[3])

# draw each little color‐tile
breaks <- seq(min(info_vals), max(info_vals), length.out=length(pal)+1)
for(i in seq_along(pal)) {
  rect(
    xleft   = x_left + (i-1) * (x_right-x_left)/length(pal),
    xright  = x_left + i     * (x_right-x_left)/length(pal),
    ybottom = y_bot,
    ytop    = y_top,
    col     = pal[i],
    border  = NA
  )
}

# add numeric labels at 5 positions
ticks <- seq(min(info_vals), max(info_vals), length.out=5)
tick_pos <- seq(x_left, x_right, length.out=5)
text(
  x    = tick_pos,
  y    = y_bot - 0.01*(usr[4]-usr[3]),
  labels = sprintf("%.2f", ticks),
  adj  = c(0.5,1),
  cex  = 0.7
)

# finally, a title for the bar
text(
  x    = (x_left + x_right)/2,
  y    = y_top + 0.01*(usr[4]-usr[3]),
  labels = "Info (bits)",
  adj  = c(0.5,0),
  cex  = 0.8
)

