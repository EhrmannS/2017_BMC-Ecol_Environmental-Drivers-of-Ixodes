################################################################################
# Introduction
#
# This script is used to conduct modelling and visualisation for an ongoing
# paper on tick abundances and how they depend on characteristics of small
# forest patches and the European agricultural landscape.
# Framework: smallFOREST
# Based on field-data collected in the framework by Steffen Ehrmann and
# field-assistants (Katja Leischke, Peter Fräßdorf, Iris Gutierrez), the
# smallFOREST site-managers and on spatial data of the smallFOREST geospatial
# database.
#
################################################################################
# License
#
# Copyright (C) 2016 Steffen Ehrmann
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful , but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# read in data and packages

setwd("/home/steffen/Documents/git/Articles/")

# read in helper functions
files <- list.files(path = "./tick_abundance_2017_smallFOREST/helper functions/", pattern = "[.]R$", recursive=T,
                    full.names = TRUE)
for(i in 1:length(files)) source(files[i])
# load in preprocessed data (the file which should be available at the location
# given in the paper)
source("./tick_abundance_2017_smallFOREST/prepare_smallFOREST-data.R")

colour1 <- "black"
colour2 <- "black"
colour3 <- "black"

################################################################################
# model building

options(contrasts=c("contr.SAS","contr.SAS"))

#---|Models|--------------------------------------------------------------------
# linear mixed model with p < 0.05 as selection criterion and R² < 0.5 as
# tolerance.
# Use: define a model, which only includes only lg.IndexAge_orig, lg.surface,
#      and the the previous tick stage. Let this run through find.effects() and
#      select the variables as described in the article.
#-------------------------------------------------------------------------------
# 1) larvae
# --->
variables <- colnames(all)[c(18:21, 28:31, 38:41, 57:dim(all)[2])]
variables <- variables[-c(which(variables=="tree_hedera_patch"),
                          which(variables=="tree_tilia_dom_patch"),
                          which(variables=="tree_medium_disp_abund"),
                          which(variables=="tree_larix_dom_patch"),
                          which(variables=="tree_pr_padus_abund"),
                          which(variables=="tree_picea_dom_patch"),
                          which(variables=="FA_structure_6"),
                          which(variables=="tree_medium_disp_alpha"),
                          which(variables=="tree_betula_dom_patch"))] # ²
l_col_lmer <- lmer(lg_l_mean ~ lg.IndexAge_orig + lg.surface + lg_a_mean + all_w_nuts_abund + I(all_w_nuts_abund^2) + temp_soil + I(temp_soil^2) + prop_past_1000 + trait13_mean + rH_5 + I(rH_5^2) + FA_structure_4 + I(FA_structure_4^2) + FA_abundances_1 + mean_alpha_tree + prop_cult_50 + herb_canopy_height_cv + ms_n.p_tot_lg_mean + I(ms_n.p_tot_lg_mean^2) + lg_shrub_weight_disp + tree_large_disp_abund + herb_regular + n.patches + I(n.patches^2) + ba_bitterlich_mean + herb_asc_to_pros + I(herb_asc_to_pros^2) +
                     (1 | id_window),
                   data = all, REML = T)
# df2 <- df; df <- find.effects(l_col_lmer,
#                               variables,
#                               degree = 2,
#                               stat = "Chisq",
#                               sbst = "p_val < 0.05 & r_sq < 0.5",
#                               order = "p_val"); View(df)
Anova(l_col_lmer, type = "III", test.statistic = "Chisq"); summary(l_col_lmer)$varcor$id_window>0
semi.residuals(l_col_lmer, c("..."), levels = "id_region")
l_col_decomp <- var.decomp(l_col_lmer, ddf="Satterthwaite", verbose = 1, penalised = T); l_col_decomp

# 2) nymphs
# --->
variables <- colnames(all)[c(5:7, 58:dim(all)[2])]
variables <- variables[-c(which(variables=="tree_hedera_patch"),
                          which(variables=="ms_p.cont_org_lg_mean"), # outlayer based trend
                          which(variables=="ms_c.cont2_lg_mean"), # outlayer based trend
                          which(variables=="tree_acer_dom_patch"), # only very little plots with value
                          which(variables=="tree_ilex_part_patch"), # only very little plots with value
                          which(variables=="tree_tilia_part_patch"), # only very little plots with value
                          which(variables=="tree_populus_part_patch"), # only little plots with value
                          which(variables=="tree_populus_dom_patch"), # only little plots with value
                          which(variables=="tree_tilia_dom_patch"))] # only very little plots with value
n_col_lmer <- lmer(lg_n_mean ~ lg.IndexAge_orig + lg.surface + lg_l_mean + ff_ph_mean + ms_c.n_lg_mean + cdd + I(cdd^2) + lg_forest.cover + prop_past_1000 + I(prop_past_1000^2) + md_total_mean + I(md_total_mean^2) + tree_w_berries_abund + rH_130 + diss_tree1_tree2 + trait11_mean + I(trait11_mean^2) + shrub_tree_total_abund + I(shrub_tree_total_abund^2) + shrub_w_nuts_abund + FA_traits_2 + lg_herb_grass + I(lg_herb_grass^2) +
                     (1 | id_window), data = all, REML = T)
df2 <- df; df <- find.effects(n_col_lmer, variables, degree = 2, stat = "Chisq", sbst = "p_val < 0.05 & r_sq < 0.5", order = "p_val"); View(df)
Anova(n_col_lmer, type = "III", test.statistic = "Chisq"); summary(n_col_lmer)$varcor$id_window>0
semi.residuals(n_col_lmer, c("..."), levels = "id_region")
# otl <- outliers(n_col_lmer, "...")
n_col_decomp <- var.decomp(n_col_lmer, ddf="Satterthwaite", verbose = 1, penalised = TRUE); n_col_decomp

# 3) adults
# --->
variables <- colnames(all)[c(5:11, 58:dim(all)[2])]
variables <- variables[-c(which(variables=="tree_hedera_patch"),
                          which(variables=="ms_p.cont_org_lg_mean"), # outlayer based trend
                          which(variables=="ms_c.cont2_lg_mean"), # outlayer based trend
                          which(variables=="tree_acer_dom_patch"), # only very little plots with value
                          which(variables=="tree_ilex_part_patch"), # only very little plots with value
                          which(variables=="tree_tilia_part_patch"), # only very little plots with value
                          which(variables=="tree_populus_part_patch"), # only little plots with value
                          which(variables=="ba_bitterlich_mean"), # only little plots with value
                          which(variables=="tree_populus_dom_patch"), # only little plots with value
                          which(variables=="FA_structure_6"), # not correlating to any other variables properly
                          which(variables=="tree_tilia_dom_patch"))] # only very little plots with value
a_col_lmer <- lmer(lg_a_mean ~ lg.IndexAge_orig + lg.surface + lg_n_mean + tree_w_nuts_e_abund + I(tree_w_nuts_e_abund^2) + shrub_w_nuts_e_abund + herb_canopy_height_cv + prop_cult_1000 + tree_total_abund + I(tree_total_abund^2) + all_large_disp_alpha + I(all_large_disp_alpha^2) + lg_prop_forest_250 + ba_bitterlich_mean + I(ba_bitterlich_mean^2) + n_med_year + I(n_med_year^2) + (1 | id_window), data = all, REML = T)
df2 <- df; df <- find.effects(a_col_lmer, variables, degree = 2, stat = "Chisq", sbst = "p_val < 0.05 & r_sq < 0.5", order = "p_val"); View(df)
Anova(a_col_lmer, type = "III", test.statistic = "Chisq"); summary(a_col_lmer)$varcor$id_window>0
semi.residuals(a_col_lmer, c("..."), levels = "id_region")
# otl <- outliers(a_col_lmer, "...")
a_col_decomp <- var.decomp(a_col_lmer, ddf="Satterthwaite", verbose = 1, penalised = TRUE); a_col_decomp


#---|response profiles|---------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) larvae
# --->
res_pro <- visreg.to.ggplot(l_col_lmer, type = "contrast", levels = "id_window", var_excl = c("id_window"))
larvae_names <- names[names(names)%in%levels(res_pro$variable)]
pos_x <- ddply(subset(res_pro, type=="fit"), .(variable), summarise, pos_x = max(x) - (max(x)-min(x))*0.5)
resp_larvae <- ggplot(subset(res_pro, type == "fit"), aes(x, y)) +
  geom_point(data = subset(res_pro, type == "res"), aes(x, y), shape = 20, size = .3) +
  geom_line(colour = colour1, size = 1) +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), bg = colour1, alpha = .3) +
  theme_bw() +
  theme(aspect.ratio=1.2,
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Log( abundance of larvae )\n") + xlab("Response variable") +
  ylim(-1.5, 2.5) +
  facet_wrap(~variable, scales = "free_x", switch = "x", drop = F, labeller = as_labeller(larvae_names), ncol = 6)

# 2) nymphs
# --->
res_pro <- visreg.to.ggplot(n_col_lmer, type = "contrast", levels = "id_window", var_excl = c("id_window"))
nymphs_names <- names[names(names)%in%levels(res_pro$variable)]
pos_x <- ddply(subset(res_pro, type=="fit"), .(variable), summarise, pos_x = max(x) - (max(x)-min(x))*0.5)
resp_nymphs <- ggplot(subset(res_pro, type == "fit"), aes(x, y)) +
  geom_point(data = subset(res_pro, type == "res"), aes(x, y), shape = 20, size = .3) +
  geom_line(colour = colour2, size = 1) +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), bg = colour2, alpha = .3) +
  theme_bw() +
  theme(aspect.ratio=1.2,
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Log( abundance of nymphs )\n") + xlab("Response variable") +
  ylim(-1, 1.5) +
  facet_wrap(~variable, scales = "free_x", switch = "x", drop = F, labeller = as_labeller(nymphs_names), ncol = 6)

# 3) adults
# --->
res_pro <- visreg.to.ggplot(a_col_lmer, type = "contrast", levels = "id_window", var_excl = c("id_window"))
adults_names <- names[names(names)%in%levels(res_pro$variable)]
pos_x <- ddply(subset(res_pro, type=="fit"), .(variable), summarise, pos_x = max(x) - (max(x)-min(x))*0.5)
resp_adults <- ggplot(subset(res_pro, type == "fit"), aes(x, y)) +
  geom_point(data = subset(res_pro, type == "res"), aes(x, y), shape = 20, size = .3) +
  geom_line(colour = colour3, size = 1) +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), bg = colour3, alpha = .3) +
  theme_bw() +
  theme(aspect.ratio=1.2,
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Log( abundance of adults )\n") + xlab("Response variable") +
  ylim(-0.9, 1.3) +
  facet_wrap(~variable, scales = "free_x", switch = "x", drop = F, labeller = as_labeller(adults_names), ncol = 6)


#---|relative importance|-------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) larvae
# --->
l_meta <- l_col_decomp$`Fixed Part`
l_meta <- cbind(l_meta, stage = "larvae")
l_meta$var2 <- as.character(rownames(l_meta))
l_meta <- join(l_meta, grouping[, c(1:2)], by = "var2")
l_meta$var1 <- ifelse(is.na(l_meta$var1), l_meta$var2, as.character(l_meta$var1))
l_meta <- join(l_meta, grouping[, c(1, 3)], by = "var1")
l_meta <- join(l_meta, grouping[, c(1, 4)], by = "var1")
l_meta <- join(l_meta, lut_group, by = "group")
l_meta <- join(l_meta, quantiles[, c(1, 3:9)], by = "var1")
colnames(l_meta)[which(colnames(l_meta)=="var1")] <- "var"
colnames(l_meta)[which(colnames(l_meta)=="var2")] <- "var_detail"
l_var <- as.data.frame(l_meta %>% group_by(var) %>% summarise(percent = round(sum(relImp_part*100, na.rm = T), 2)))
l_var$label <- sprintf("eta**2==%.1f", l_var$percent)
colnames(l_var)[which(colnames(l_var)=="var")] <- "variable"
l_var <- join(l_var, pos_x, by = "variable")
l_var$variable <- factor(l_var$variable, levels(res_pro$variable)); l_var <- l_var[order(l_var$variable),]
l_drivergroup <- l_meta %>%
  group_by(metagroup) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); l_drivergroup
l_group <- l_meta %>%
  group_by(group) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); l_group
l_habitat <- l_meta %>%
  filter(type == "microhabitat" | type == "macrohabitat") %>%
  group_by(type) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); l_habitat
resp_larvae <- resp_larvae + geom_text(data = l_var, aes(x = pos_x, y = 2, label = label), parse = T)
svg(filename = "resp_larvae.svg", height = 20, width = 11); resp_larvae; dev.off()

# 2) nymphs
# --->
n_meta <- n_col_decomp$`Fixed Part`
n_meta <- cbind(n_meta, stage = "nymph")
n_meta$var2 <- as.character(rownames(n_meta))
n_meta <- join(n_meta, grouping[, c(1:2)], by = "var2")
n_meta$var1 <- ifelse(is.na(n_meta$var1), n_meta$var2, as.character(n_meta$var1))
n_meta <- join(n_meta, grouping[, c(1, 3)], by = "var1")
n_meta <- join(n_meta, grouping[, c(1, 4)], by = "var1")
n_meta <- join(n_meta, lut_group, by = "group")
n_meta <- join(n_meta, quantiles[, c(1, 3:9)], by = "var1")
colnames(n_meta)[which(colnames(n_meta)=="var1")] <- "var"
colnames(n_meta)[which(colnames(n_meta)=="var2")] <- "var_detail"
n_var <- as.data.frame(n_meta %>% group_by(var) %>% summarise(percent = round(sum(relImp_part*100, na.rm = T), 2)))
n_var$label <- sprintf("eta**2==%.1f", n_var$percent)
colnames(n_var)[which(colnames(n_var)=="var")] <- "variable"
n_var <- join(n_var, pos_x, by = "variable")
n_var$variable <- factor(n_var$variable, levels(res_pro$variable)); n_var <- n_var[order(n_var$variable),]
n_drivergroup <- n_meta %>%
  group_by(metagroup) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); n_drivergroup
n_group <- n_meta %>%
  group_by(group) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); n_group
n_habitat <- n_meta %>%
  filter(type == "microhabitat" | type == "macrohabitat") %>%
  group_by(type) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); n_habitat
resp_nymphs <- resp_nymphs + geom_text(data = n_var, aes(x = pos_x, y = 1.2, label = label), parse = T)
svg(filename = "resp_nymphs.svg", height = 20, width = 11); resp_nymphs; dev.off()

# 3) adults
# --->
a_meta <- a_col_decomp$`Fixed Part`
a_meta <- cbind(a_meta, stage = "adult")
a_meta$var2 <- as.character(rownames(a_meta))
a_meta <- join(a_meta, grouping[, c(1:2)], by = "var2")
a_meta$var1 <- ifelse(is.na(a_meta$var1), a_meta$var2, as.character(a_meta$var1))
a_meta <- join(a_meta, grouping[, c(1, 3)], by = "var1")
a_meta <- join(a_meta, grouping[, c(1, 4)], by = "var1")
a_meta <- join(a_meta, lut_group, by = "group")
a_meta <- join(a_meta, quantiles[, c(1, 3:9)], by = "var1")
colnames(a_meta)[which(colnames(a_meta)=="var1")] <- "var"
colnames(a_meta)[which(colnames(a_meta)=="var2")] <- "var_detail"
a_var <- as.data.frame(a_meta %>% group_by(var) %>% summarise(percent = round(sum(relImp_part*100, na.rm = T), 2)))
a_var$label <- sprintf("eta**2==%.1f", a_var$percent)
colnames(a_var)[which(colnames(a_var)=="var")] <- "variable"
a_var <- join(a_var, pos_x, by = "variable")
a_var$variable <- factor(a_var$variable, levels(res_pro$variable)); a_var <- a_var[order(a_var$variable),]
a_drivergroup <- a_meta %>%
  group_by(metagroup) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); a_drivergroup
a_group <- a_meta %>%
  group_by(group) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); a_group
a_habitat <- a_meta %>%
  filter(type == "microhabitat" | type == "macrohabitat") %>%
  group_by(type) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); a_habitat
resp_adults <- resp_adults + geom_text(data = a_var, aes(x = pos_x, y = 1, label = label), parse = T)
svg(filename = "resp_adults.svg", height = 20, width = 11); resp_adults; dev.off()


################################################################################
# Correlations

correlations <- ddply(all, .(id_window), summarise,
                      Y_lat = mean(Y_lat),
                      doy = mean(doy),
                      gdd = mean(gdd),
                      cdd = mean(cdd),
                      precip_pre = mean(precip_pre),
                      precip_year = mean(precip_year),
                      temp_mean_pre = mean(temp_mean_pre),
                      temp_mean_year = mean(temp_mean_year),
                      temp_min_pre = mean(temp_min_pre),
                      temp_min_year = mean(temp_min_year),
                      temp_max_pre = mean(temp_max_pre),
                      temp_max_year = mean(temp_max_year),
                      n_snow_year = mean(n_snow_year),
                      n_dew_pre = mean(n_dew_pre),
                      n_dew_year = mean(n_dew_year),
                      n_cold_pre = mean(n_cold_pre),
                      n_cold_year = mean(n_cold_year),
                      n_med_pre = mean(n_med_pre),
                      n_med_year = mean(n_med_year),
                      n_warm_pre = mean(n_warm_pre),
                      n_warm_year = mean(n_warm_year)
)

datas <- as.matrix(correlations[,c(2:dim(correlations)[2])])
corrs <- rcorr(datas)
plotcorr(corr = corrs, outline = TRUE)

################################################################################
# Graphs

tick_stats <- rbind(l_meta, n_meta, a_meta)
tick_stats$group <- droplevels(tick_stats$group)
tick_stats$group[which(tick_stats$group=="taxonomic properties")] <- "functional properties"
write.csv(tick_stats, "./tick_abundance_2017_smallFOREST/tick-models.csv")

# relative importance of variable groups
tick_group <- ddply(tick_stats, .(stage, group), summarise,
                    relImp_part = sum(relImp_part*100))

tick_group <- reshape(tick_group, v.names = "relImp_part", timevar = "group", idvar = c("stage"), direction = "wide")
tick_group[is.na(tick_group)] <- 0.01
tick_group <- reshape(tick_group, direction = "long"); row.names(tick_group) <- NULL

tick_group$group <- factor(tick_group$group, levels(tick_group$group)[c(3, 2, 1, 6, 4, 8, 5)]); tick_group <- tick_group[order(tick_group$stage, tick_group$group),]

# relative importance of variable meta groups
tick_metagroup <- ddply(tick_stats, .(stage, metagroup), summarise,
                        relImp_part = sum(relImp_part*100))
tick_metagroup <- reshape(tick_metagroup, v.names = "relImp_part", timevar = "metagroup", idvar = c("stage"), direction = "wide")
tick_metagroup[is.na(tick_metagroup)] <- 0
tick_metagroup <- reshape(tick_metagroup, direction = "long"); row.names(tick_metagroup) <- NULL

tick_metagroup$metagroup <- factor(tick_metagroup$metagroup, levels(tick_metagroup$metagroup)[c(3, 2, 1, 4)]); tick_metagroup <- tick_metagroup[order(tick_metagroup$stage, tick_metagroup$metagroup),]

tick_metagroup_avrg <- ddply(tick_metagroup, .(metagroup), summarise,
                         relImp_avrg = mean(relImp_part))

g1 <- ggplot(tick_metagroup, aes(metagroup, relImp_part, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge", size = .3) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("\nDriver group") + ylab("Relative importance [%]\n") +
  expand_limits(y=c(0,70)) +
  scale_fill_manual(name = "Tick Stage", values = c("gray90", "gray60", "gray30"), labels = c("Larvae", "Nymphs", "Adults")) +
  scale_x_discrete(labels = c("Macroclimate\n", "Landscape", "Habitat", "Ontogeny")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

g1_avrg <- ggplot(tick_metagroup_avrg, aes(metagroup, relImp_avrg)) +
  geom_bar(stat = "identity", size = .3) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("\nDriver group") + ylab("Average relative importance [%]\n") +
  expand_limits(y=c(0,70)) +
  scale_x_discrete(labels = c("Macroclimate\n", "Landscape", "Habitat", "Ontogeny")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

# relative importance of only habitat variables
tick_stats_sub <- subset(tick_stats, subset = tick_stats$metagroup=="habitat")
tick_stats_sub <- droplevels(tick_stats_sub)

tick_group2 <- ddply(tick_stats_sub, .(stage, group), summarise,
                     relImp_part = sum(relImp_part*100))
tick_group3 <- ddply(tick_stats_sub, .(stage, type), summarise,
                     relImp_part = sum(relImp_part*100))


tick_group2 <- reshape(tick_group2, v.names = "relImp_part", timevar = "group", idvar = c("stage"), direction = "wide")
tick_group2[is.na(tick_group2)] <- 0.01
tick_group2 <- reshape(tick_group2, direction = "long"); row.names(tick_group2) <- NULL

tick_group2$group <- factor(tick_group2$group, levels(tick_group2$group)[c(1, 4, 2, 3)]); tick_group2 <- tick_group2[order(tick_group2$stage, tick_group2$group),]
tick_group2_avrg <- ddply(tick_group2, .(group), summarise,
                          relImp_avrg = mean(relImp_part))

tick_group3$type <- factor(tick_group3$type, levels(tick_group3$type)[c(2, 1)]); tick_group3 <- tick_group3[order(tick_group3$stage, tick_group3$type),]
tick_group3_avrg <- ddply(tick_group3, .(type), summarise,
                          relImp_avrg = mean(relImp_part))

g2 <- ggplot(tick_group3, aes(type, relImp_part, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge", size = .3) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("\n...within Habitat") +
  expand_limits(y=c(0,70)) +
  scale_fill_manual(name = "Tick Stage", values = c("gray90", "gray60", "gray30"), labels = c("Larvae", "Nymphs", "Adults")) +
  scale_x_discrete(labels = c("Microhabitat\n", "Macrohabitat")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

g2_avrg <- ggplot(tick_group3_avrg, aes(type, relImp_avrg)) +
  geom_bar(stat = "identity", size = .3) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("\n...within Habitat") +
  expand_limits(y=c(0,70)) +
  scale_x_discrete(labels = c("Microhabitat\n", "Macrohabitat\n")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

g3 <- ggplot(tick_group2, aes(group, relImp_part, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge", size = .3) +
  theme_bw() +
  xlab("\n...within Habitat") +
  expand_limits(y=c(0,70)) +
  scale_fill_manual(name = "Tick Stage", values = c("gray90", "gray60", "gray30"), labels = c("Larvae", "Nymphs", "Adults")) +
  scale_x_discrete(labels = c("Functional\nproperties", "Structural\nproperties", "Microclimate", "Soil")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.background = element_rect(colour = "lightgray"),
        legend.position = c(.85,.85),
        legend.key = element_rect(colour = "white"))

g3_avrg <- ggplot(tick_group2_avrg, aes(group, relImp_avrg)) +
  geom_bar(stat = "identity", size = .3) +
  theme_bw() +
  xlab("\n...within Habitat") +
  expand_limits(y=c(0,70)) +
  scale_x_discrete(labels = c("Functional\nproperties", "Structural\nproperties", "Microclimate", "Soil")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.background = element_rect(colour = "lightgray"),
        legend.position = c(.85,.85),
        legend.key = element_rect(colour = "white"))

multiplot(g1, g2, g3, cols = 3, layout = matrix(c(1, 1, 1, 1, 2, 2, 3, 3, 3, 3), nrow = 1))
multiplot(g1_avrg, g2_avrg, g3_avrg, cols = 3, layout = matrix(c(1, 1, 1, 1, 2, 2, 3, 3, 3, 3), nrow = 1))

################################################################################
# footnotes

# ¹ exclude outlaying values, rerun model and decide based on this to include
# "significant" variable, or not...
# ² look at semiresidual plot to find strange patterns and significances based
# on outlayers.
