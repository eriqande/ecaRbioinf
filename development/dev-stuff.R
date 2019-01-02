
combo <- haps$tidy

smat <- combo %>%
  unite(col = "name", sep = "-", Indiv, haplo) %>%
  select(-ChromKey) %>%
  spread(key = name, value = allele) %>%
  select(-POS) %>%
  as.matrix() %>%
  t()


D <- phyDat(smat)


d <- dist.ml(D)
nj <- NJ(d)
plot(nj, )
write.tree(nj, file = "~/Desktop/tree.txt")

######################################################

rosa <- haps$avd %>% filter(POS > 12.05e6, POS < 12.4e06)

boing <- ggplot(rosa, aes(x = factor(POS), y = Indiv, fill = anc_vs_derived)) +
  geom_raster() +
  scale_fill_manual(values = c(D = "yellow", A = "blue"))

ggsave(boing, filename = "/tmp/raster.jpg", width = 32, height = 20)
