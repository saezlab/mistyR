library(tidyverse)
library(cowplot)
library(scales)
library(igraph)


# make a shiny app? control diferent values

# images <- list.files(path = "melanoma_results/ROI", "*PRE", full.names = TRUE)[-6]
images <- list.dirs(path = "bc_results/", full.names = TRUE)[-1]
# images <- list.files(path = "synthetic_results/", full.names = TRUE)

# global p-value significance cutoff
p.cutoff <- 0.05


impr <- images %>% map_dfc(function(image) {
  performance <- read_delim(paste0(image, .Platform$file.sep, "performance.txt"),
    delim = " ", col_types = cols()
  ) %>% distinct()

  targets <- unique(performance$target)
  performance %>%
    arrange(target) %>%
    transmute(RMSE = (intra.RMSE - multi.RMSE) / intra.RMSE, R2 = (multi.R2 - intra.R2))
})

# RMSE
ggplot(impr %>% select(contains("RMSE")) %>%
  mutate(target = targets) %>%
  pivot_longer(cols = -target, names_to = "name", values_to = "value")) +
  geom_boxplot(aes(x = target, y = value * 100)) +
  theme_classic() +
  ylab("Improvement (%)") +
  xlab("Target") +
  ylim(c(-5, 17)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

insig.index.RMSE <- which(impr %>%
  select(contains("RMSE")) %>%
  apply(1, t.test, alternative = "greater") %>%
  map_dbl(~ .x$p.value) %>% p.adjust(method = "fdr") > p.cutoff)

# R2
ggplot(impr %>% select(contains("R2")) %>%
  mutate(target = targets) %>%
  pivot_longer(cols = -target, names_to = "name", values_to = "value") %>%
  filter(value > -1)) +
  geom_boxplot(aes(x = target, y = value * 100)) +
  theme_classic() +
  ylab("Change in variance explained") +
  xlab("Target") +
  ylim(c(-20, 25)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

insig.index.R2 <- which(impr %>%
  select(contains("R2")) %>%
  apply(1, t.test, alternative = "greater") %>%
  map_dbl(~ .x$p.value) %>%
  p.adjust(method = "fdr") > p.cutoff)

avg <- ((images %>% map(function(image) {
  coefficients <- read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
    delim = " ", col_types = cols()
  )

  targets <<- coefficients %>%
    pull(target) %>%
    sort

  coefficients %>%
    distinct() %>%
    arrange(target) %>%
    select(-target, -contains("intercept")) %>%
    mutate_at(vars(starts_with("p.")), ~ as.numeric(. <= p.cutoff)) %>%
    mutate_at(vars(-starts_with("p.")), abs)
}) %>% 
  reduce(`+`)) / length(images)) %>% 
  mutate(target = targets)

ctotals <- avg %>% 
  select(-starts_with("p."), -"target") %>%
  rowSums

coefs <- avg %>% 
  select(-starts_with("p.")) %>%
  mutate_if(is.numeric, ~./ctotals) %>%
  pivot_longer(-target, names_to = "view")

ggplot(coefs) + 
  geom_col(aes(x=target, y=value, group=view, fill=view)) +
  scale_fill_brewer(palette="Dark2") +
  xlab("Target") +
  ylab("Contribution") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

maps <- images %>% map(function(image) {
  coefficients <- read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
    delim = " ", col_types = cols()
  ) %>% distinct()

  targets <- unique(coefficients$target)
  views <<- (coefficients %>% select(-target, -starts_with("p."), -intercept) %>% colnames())

  # one heatmap per view
  maps <- views %>% map(function(view) {
    all.importances <- targets %>% map(~ read_csv(paste0(
      image, .Platform$file.sep, "importances_",
      .x, "_", view, ".txt"
    ),
    col_types = cols()
    ) %>%
      distinct() %>%
      filter(!grepl("_2$", target)))

    features <- unique(all.importances %>% map(~ .x$target) %>% unlist())

    pview <- paste0("p.", view)
    ps <- coefficients %>%
      select(target, !!pview) %>%
      mutate(!!pview := (1 - !!sym(pview)))


    # importances are standardized for each target an multiplied by 1-pval(view)
    result <- all.importances %>%
      imap_dfc(~
      tibble(target = features, zero.imp = 0) %>%
        left_join(.x, by = "target") %>%
        transmute(feature = target, importance = (zero.imp + scale(imp)[, 1]) *
          (ps %>% filter(target == targets[.y]) %>% pull(pview))) %>%
        select(importance)) %>%
      `colnames<-`(targets) %>%
      mutate(Predictor = features)

    # in order for aggregation
    result %>%
      arrange(Predictor) %>%
      select(noquote(order(colnames(result))))
  })
})

aggregated <- maps %>% reduce(function(acc, l) {
  map2(acc, l, ~ (((.x %>% select(-Predictor)) + (.y %>% select(-Predictor))) %>%
    mutate(Predictor = .x %>% pull(Predictor))))
})

# set all rows to zero where q > cutoff?

# aggregated <- aggregated %>% map(function(view){
#   view[insig.index.RMSE, -27] <- 0
#   view
# })



heatmaps <- aggregated %>% imap(function(hmap, view.index) {
  ggplot(hmap %>% gather(key = "Target", value = "importance", -Predictor)) +
    geom_tile(aes(x = Predictor, y = Target, fill = importance / length(images))) +
    scale_fill_gradient2(low = "white", mid = "white", high = muted("blue"), midpoint = 0) +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle(views[view.index])
})

plot_grid(plotlist = heatmaps, labels = "AUTO", align = "hv")


# importance threshold should be a parameter
imp.thresh <- 0.5

par(mfrow = c(2, 2), mai = c(0.1, 0.1, 0.5, 0.1))
communities <- aggregated %>% imap(function(m, i) {
  A <- (m %>% arrange(Predictor) %>% select(sort(m %>% pull(Predictor))) %>% as.matrix()) / length(images)
  # filter negative values and make the diagonal 0
  A[A < imp.thresh | is.na(A)] <- 0

  # custom filter for melanoma study names
  G <- graph.adjacency(A, mode = "plus", weighted = TRUE) %>%
    set.vertex.attribute("name", value = gsub("Cell_[0-9]+(PRE|ON|POST)", "", names(V(.))))

  C1 <- cluster_louvain(G)

  layout <- layout.fruchterman.reingold(G)
  plot(G, layout = layout, mark.groups = C1, main = views[i], vertex.size = 4, vertex.color = "black", vertex.label.dist = 1)
})
