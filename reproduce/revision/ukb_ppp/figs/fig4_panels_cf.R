## Revision copy of code/four_panels_and_cis_trans.R (Fig 4 panels and hotspot plots) reading the
## no-PRS / no-PC cross-fit results (results/web_SynGWAS_summary_results_cf_noPRS_noPC/) and
## writing to results/figs_cf_noPRS_noPC/.
setwd("<PROJECT_DIR>")
library(readr)
#protein_info <- read_csv("data/protein_info.csv")
#View(protein_info)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

#### 1. find the panel and summary ##########
final_protein_map <- readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds")
unique(final_protein_map$panel)

# Directory and filename pattern
dir_path1 <- "<PROJECT_DIR>/results/table_results"
dir_path2 <- "<PROJECT_DIR>/results/web_SynGWAS_summary_results_cf_noPRS_noPC"

file_path <- file.path(dir_path2, "significant_loci_v2geno_all_crossfit_noPRS.csv")
df <- read.csv(file_path)

num_cols <- c("threshold", "predicted_correlation",
              "SynGWAS_hat_beta", "SynGWAS_sd_hat_beta", "SynGWAS_pvalues",
              "ObsGWAS_hat_beta", "ObsGWAS_sd_hat_beta", "ObsGWAS_pvalues")

df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

thresh_vals <- sort(unique(df$threshold),decreasing = T)

counts <- df %>%
  group_by(protein_short_name, predicted_correlation) %>%
  summarise(
    threshold  = list(thresh_vals),
    syn_count  = list(sapply(thresh_vals, function(t) sum(SynGWAS_pvalues  <= t, na.rm = TRUE))),
    obs_count  = list(sapply(thresh_vals, function(t) sum(ObsGWAS_pvalues  <= t, na.rm = TRUE))),
    .groups = "drop"
  ) %>% unnest(c(threshold, syn_count, obs_count))

counts_5e08 <- counts[counts$threshold==thresh_vals[1],]
## revision: the counts above add up the rows of all four thresholds and the union of
## Syn-PALM and Obs-GWAS lead variants, so a locus was counted up to four times; use
## the numbers of primary associations of each method at 5e-8 instead (aggregate_all.R)
.per <- read.csv("<PROJECT_DIR>/results/table_results_cf_noPRS_noPC/summary_per_protein.csv", check.names = FALSE)
counts_5e08 <- data.frame(protein_short_name = .per$protein, predicted_correlation = .per$rho_oof, threshold = 5e-8,
                          syn_count = .per[["SynG_loci_5e-08"]], obs_count = .per[["ObsG_loci_5e-08"]])
counts_5e08$sign_diff <- counts_5e08$syn_count-counts_5e08$obs_count
# View(counts_5e08)


counts_5e08$protein_short_name <- gsub("-", "_", counts_5e08$protein_short_name)
all(counts_5e08$protein_short_name %in% final_protein_map$id)

counts_with_panel <- counts_5e08 %>% left_join(final_protein_map %>% select(id, panel), by = c("protein_short_name" = "id"))

counts_with_panel$panel4 <- counts_with_panel$panel
library(stringr)
library(dplyr)

counts_with_panel <- counts_with_panel %>% 
  mutate(panel4 = str_remove(panel4, "_II"))

head(counts_with_panel)


library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

# --- 1. 聚合数据并计算差值标签 ---
plot_data_panel <- counts_with_panel %>%
  group_by(panel4) %>%
  summarise(
    Synthetic = sum(syn_count, na.rm = TRUE),
    Observed = sum(obs_count, na.rm = TRUE),
    n_prot = n(),
    .groups = 'drop'
  ) %>%
  # revision (Reviewer 3, minor 3): label with the number of proteins and the
  # additional associations per protein instead of the total difference
  mutate(diff = Synthetic - Observed) %>%
  mutate(x_label = sprintf("%s\n(%s proteins; +%.2f per protein)", panel4, scales::comma(n_prot), diff / n_prot))

# --- 2. 颜色定义 ---
my_colors_panel <- c(
  "Inflammation"    = "#E79D4F", 
  "Cardiometabolic" = "#31A231", 
  "Oncology"        = "#3274B5", 
  "Neurology"       = "#8B6BB2"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

# --- 1. 准备并列绘图的长格式数据 ---
plot_data_long <- plot_data_panel %>%
  # 将 Synthetic 和 Observed 转为一列，方便 position_dodge
  pivot_longer(cols = c(Synthetic, Observed), 
               names_to = "Method_Type", 
               values_to = "Count_Value") %>%
  # 统一图例名称
  mutate(Method_Type = case_when(
    Method_Type == "Synthetic" ~ "Syn-PALM",
    Method_Type == "Observed" ~ "Obs-GWAS"
  ))

# --- 绘图 ---
panel4 <- ggplot(plot_data_long, aes(x = fct_reorder(x_label, Count_Value, .fun = max), 
                           y = Count_Value, 
                           fill = panel4, 
                           alpha = Method_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = comma(Count_Value)), 
            position = position_dodge(width = 0.8), 
            hjust = 1.2,              # 关键：正值会让文字向左（柱子内部）移动
            size = 3, 
            fontface = "bold", 
            color = "black",          # 内部显示通常白色更清晰
            show.legend = FALSE) +
  scale_fill_manual(values = my_colors_panel) +
  scale_alpha_manual(
    name = "Method:", 
    values = c("Obs-GWAS" = 0.3, "Syn-PALM" = 1), 
    breaks = c("Obs-GWAS", "Syn-PALM")
  ) +
  coord_flip() +
  # 内部显示不需要太大的右侧 expand
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), labels = comma) +
  
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = "black", size = 10, face = "bold"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  ) +
  
  # 5. 隐藏 fill 图例
  guides(
    fill = "none",
    alpha = guide_legend(
      override.aes = list(
        fill = "grey50", 
        alpha = c(0.3, 1)
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Sum of GWS Associations"
  )

ggsave(
  plot = panel4,
  file = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/panel4_results.pdf",
  width = 4.9,    # revision 2026-10-19: smaller canvas for the compact Fig 4 layout (was 6.0 x 3.2)
  height = 2.6,
  units = "in",
  dpi = 480
)


#### 2. find the disease and summary ##########
final_protein_map <- readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds")

# Directory and filename pattern
dir_path1 <- "<PROJECT_DIR>/results/table_results"
dir_path2 <- "<PROJECT_DIR>/results/web_SynGWAS_summary_results_cf_noPRS_noPC"

file_path <- file.path(dir_path2, "significant_loci_v2geno_all_crossfit_noPRS.csv")
df <- read.csv(file_path)

num_cols <- c("threshold", "predicted_correlation",
              "SynGWAS_hat_beta", "SynGWAS_sd_hat_beta", "SynGWAS_pvalues",
              "ObsGWAS_hat_beta", "ObsGWAS_sd_hat_beta", "ObsGWAS_pvalues")

df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

thresh_vals <- sort(unique(df$threshold),decreasing = T)

counts <- df %>%
  group_by(protein_short_name, predicted_correlation) %>%
  summarise(
    threshold  = list(thresh_vals),
    syn_count  = list(sapply(thresh_vals, function(t) sum(SynGWAS_pvalues  <= t, na.rm = TRUE))),
    obs_count  = list(sapply(thresh_vals, function(t) sum(ObsGWAS_pvalues  <= t, na.rm = TRUE))),
    .groups = "drop"
  ) %>% unnest(c(threshold, syn_count, obs_count))

counts_5e08 <- counts[counts$threshold==thresh_vals[1],]
## revision: the counts above add up the rows of all four thresholds and the union of
## Syn-PALM and Obs-GWAS lead variants, so a locus was counted up to four times; use
## the numbers of primary associations of each method at 5e-8 instead (aggregate_all.R)
.per <- read.csv("<PROJECT_DIR>/results/table_results_cf_noPRS_noPC/summary_per_protein.csv", check.names = FALSE)
counts_5e08 <- data.frame(protein_short_name = .per$protein, predicted_correlation = .per$rho_oof, threshold = 5e-8,
                          syn_count = .per[["SynG_loci_5e-08"]], obs_count = .per[["ObsG_loci_5e-08"]])
counts_5e08$sign_diff <- counts_5e08$syn_count-counts_5e08$obs_count
# View(counts_5e08)

counts_5e08$protein_short_name <- gsub("-", "_", counts_5e08$protein_short_name)
all(counts_5e08$protein_short_name %in% final_protein_map$id)

###### load diease info
protein_disease <- read_csv("data/protein_disease.csv")

main_headers <- colnames(protein_disease)
sub_headers <- as.character(protein_disease[1, ])
current_main <- ""
for(i in seq_along(main_headers)) {
  if(!is.na(main_headers[i]) && !str_detect(main_headers[i], "^\\.\\.\\.")) {
    current_main <- main_headers[i]
  }
  main_headers[i] <- current_main
}

new_colnames <- paste(main_headers, sub_headers, sep = "_")

new_colnames[1] <- "UKBPPP_ProteinID"
new_colnames[2] <- "Protein_name"

colnames(protein_disease) <- new_colnames
protein_disease_clean <- protein_disease[-1, ] # 去掉原本的第一行

# 6. 按照你之前的要求：提取 ProteinID 并转小写
protein_disease_clean <- protein_disease_clean %>%
  filter(!is.na(UKBPPP_ProteinID)) %>%
  mutate(protein_short_name = tolower(sub(":.*", "", UKBPPP_ProteinID)))

# 检查结果
head(protein_disease_clean[, 1:10])
protein_disease_clean$protein_short_name <- gsub("-", "_", protein_disease_clean$protein_short_name)
all(protein_disease_clean$protein_short_name %in% final_protein_map$id)

disease_used <- unique(main_headers)[-c(1:2,28)]

thre <- 15
protein_disease_num <- protein_disease_clean %>% mutate(across(contains("log10(p)"), as.numeric))

library(dplyr)
library(purrr)

# 1. 将 counts_5e08 的 protein_short_name 统一转为小写（确保匹配成功）
counts_5e08 <- counts_5e08 %>%
  mutate(protein_short_name = tolower(protein_short_name))

library(dplyr)
library(stringr)

# 存储完整结果的列表
results_list <- list()
# 专门存储 Top1 蛋白的名字和对应的 P 值（用于后续画图标签）
top_protein_df <- data.frame()

for (disease in disease_used) {
  
  # 动态匹配当前疾病的 log10(p) 列名
  target_col <- colnames(protein_disease_num) %>% 
    str_subset(fixed(disease)) %>% 
    str_subset("log10\\(p\\)")
  
  if (length(target_col) > 0) {
    # 提取该疾病下所有显著的行 (log10P > 10)
    sig_data <- protein_disease_num %>%
      filter(!!sym(target_col) > thre) %>%
      arrange(desc(!!sym(target_col))) # 按显著性降序排列
    
    # 存储到 results_list
    results_list[[disease]] <- sig_data$protein_short_name
    
    # 如果有显著蛋白，记录第一名
    if (nrow(sig_data) > 0) {
      top_info <- data.frame(
        disease = disease,
        top_protein = sig_data$protein_short_name[1],
        top_log10P = sig_data[[target_col]][1],
        stringsAsFactors = FALSE
      )
      top_protein_df <- rbind(top_protein_df, top_info)
    }
  }
}
print(results_list[["Smoking"]])

disease_counts <- sapply(results_list, length)
print(sort(disease_counts, decreasing = TRUE))

# 2. 使用 map_df 遍历列表并生成汇总行
disease_summary <- map_df(names(results_list), function(disease_name) {
  
  # 获取当前疾病包含的蛋白名字列表
  proteins <- results_list[[disease_name]]
  
  # 在 counts_5e08 中匹配这些蛋白
  matched_data <- counts_5e08 %>%
    filter(protein_short_name %in% proteins)
  
  # 计算各项指标
  data.frame(
    disease = disease_name,
    protein_count = length(proteins),             # 原始列表中的蛋白数
    matched_protein_count = nrow(matched_data),    # 在 counts 表中实际匹配到的数
    total_syn_count = sum(matched_data$syn_count, na.rm = TRUE),
    total_obs_count = sum(matched_data$obs_count, na.rm = TRUE)
  )
})

# 3. 查看汇总表
print(disease_summary)

## Reviewer 3 minor 3: category sizes, associations per protein and enrichment of
## the additional Syn-PALM associations (Fisher test of the group against all
## other proteins: additional vs Obs-GWAS associations)
cnt <- as.data.frame(counts_5e08)
cnt$extra <- pmax(cnt$syn_count - cnt$obs_count, 0)
disease_rates <- map_df(names(results_list), function(d) {
  m <- cnt[cnt$protein_short_name %in% results_list[[d]], ]
  r <- cnt[!cnt$protein_short_name %in% results_list[[d]], ]
  data.frame(disease = d, proteins = nrow(m), syn = sum(m$syn_count), obs = sum(m$obs_count),
             syn_per_protein = mean(m$syn_count), obs_per_protein = mean(m$obs_count),
             extra_per_protein = mean(m$extra), ratio = sum(m$syn_count) / sum(m$obs_count),
             fisher_p = fisher.test(matrix(c(sum(m$extra), sum(m$obs_count), sum(r$extra), sum(r$obs_count)), 2))$p.value,
             odds_ratio = unname(fisher.test(matrix(c(sum(m$extra), sum(m$obs_count), sum(r$extra), sum(r$obs_count)), 2))$estimate))
})
disease_rates <- disease_rates[order(-disease_rates$extra_per_protein), ]
write.csv(disease_rates, "<PROJECT_DIR>/results/table_results_cf_noPRS_noPC/disease_group_rates.csv", row.names = FALSE)
print(disease_rates)



### 画一部分disease吧
disease_summary_sub <- disease_summary[order(disease_summary$total_syn_count - disease_summary$total_obs_count,decreasing = T)[1:15],]

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# --- Step 1: 准备最全的自定义颜色表 (Named Vector) ---
# 确保 disease_used 中的每一个字符串在这里都有对应的颜色
my_disease_colors <- c(
  "ALT (liver function)" = "#984EA3",                
  "AST (liver function)" = "#4DAF4A",           
  "eGFR (renal function)" = "#E41A1C",
  "Medication use" = "#56B4E9",
  "Smoking" = "#E69F00", 
  "Essential (primary) hypertension" = "#CC79A7",
  "Chronic ischaemic heart disease" = "#0072B2",
  "Disorders of lipoprotein metabolism and other lipidaemias" = "#009E73",
  "Gastritis and duodenitis" = "#FFD92F",
  "Unspecified acute lower respiratory infection" = "#CCEEFF",
  "Other diseases of anus and rectum" = "#FC8D62",
  "Depressive episode" = "#C009A7",
  "Asthma" = "#A6D854",
  "Diaphragmatic hernia" = "#66C2A5",
  "Diverticular disease of intestine" = "#8DA0CB"
)

library(ggplot2)
library(scales)
library(dplyr)
library(stringr)

plot_summary_final <- disease_summary_sub %>%
  # 1. 计算差值：Synthetic - Observed
  mutate(diff_count = total_syn_count - total_obs_count) %>%
  # 2. 构造新的 X 轴标签：疾病名 (差值)
  # 这里使用了 comma 函数，确保如果数字很大，会有千分位分隔符（如 1,234）
  mutate(x_label = sprintf("%s (%s proteins; +%.2f per protein)", disease, scales::comma(matched_protein_count), diff_count / matched_protein_count)) %>%
  left_join(top_protein_df, by = "disease") %>%
  arrange(total_syn_count)

plot_summary_long <- plot_summary_final %>%
  select(
    disease,
    x_label,
    top_protein,
    total_obs_count,
    total_syn_count
  ) %>%
  pivot_longer(
    cols = c(total_obs_count, total_syn_count),
    names_to = "Method_Type",
    values_to = "pQTL_count"
  ) %>%
  mutate(
    Method_Type = recode(
      Method_Type,
      total_obs_count = "Obs-GWAS",
      total_syn_count = "Syn-PALM"
    ),
    x_label_wrapped = str_wrap(x_label, width = 40)
  )

# --- Step 3: 绘图 ---
library(ggplot2)
library(scales)
library(dplyr)
library(stringr)

library(stringr)

# --- 修改数据预处理部分 ---
plot_summary_long <- plot_summary_long %>%
  # 对 x_label 进行换行处理，width = 30 表示大约30个字符换一行
  # 这里的 reorder 逻辑要保持一致
  mutate(x_label_wrapped = str_wrap(x_label, width = 40))

# --- 修改绘图代码中的映射 ---
plot_disease10 <- ggplot(plot_summary_long, 
                         aes(x = reorder(x_label_wrapped, pQTL_count, .fun = max), # 使用换行后的标签
                             y = pQTL_count, 
                             fill = disease, 
                             alpha = Method_Type)) +
  
  # 1. 绘制并列柱子
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  # 2. 数字在 bar 内部
  geom_text(aes(label = comma(pQTL_count)), 
            position = position_dodge(width = 0.8), 
            hjust = 1.2,              
            size = 1.1,                
            fontface = "bold", 
            color = "black",           
            show.legend = FALSE) +
  
  # 3. Top Protein 在 bar 外部
  geom_text(data = subset(plot_summary_long, Method_Type == "Syn-PALM"),
            aes(label = top_protein), 
            position = position_dodge(width = 0.8), 
            hjust = -0.15,             
            size = 2.8, 
            fontface = "italic", 
            color = "black",
            show.legend = FALSE) +
  
  # 4. 颜色与图例配置
  scale_fill_manual(values = my_disease_colors) + 
  scale_alpha_manual(
    name = "Association Method:", 
    values = c("Obs-GWAS" = 0.3, "Syn-PALM" = 1), 
    breaks = c("Obs-GWAS", "Syn-PALM")
  ) +
  
  # 5. 坐标轴翻转与比例调整
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13)), labels = comma) + 
  
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    # 关键：设置 lineheight（行间距），避免换行后的文字挤在一起
    axis.text.y = element_text(color = "black", size = 7, lineheight = 0.9), 
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 8), 
    legend.key = element_blank(),
    legend.margin = margin(t = -5)
  ) +
  
  # 6. 图例美化
  guides(
    fill = "none", 
    alpha = guide_legend(
      override.aes = list(
        fill = "grey50", 
        alpha = c(0.3, 1)
      )
    )
  ) + 
  labs(
    x = NULL,
    y = "Number of significant associations"
  )

print(plot_disease10)

ggsave(
  plot = plot_disease10,
  file = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/plot_disease10_results.pdf",
  width = 7.0,
  height = 4.5,
  units = "in",
  dpi = 680
)


#### 3. from the pos to summary cis and trans  ##########

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ==========================================
# 1. 导入与基础清洗
# ==========================================
final_protein_map <- readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds")

dir_path2 <- "<PROJECT_DIR>/results/web_SynGWAS_summary_results_cf_noPRS_noPC"
file_path <- file.path(dir_path2, "significant_loci_v2geno_all_crossfit_noPRS.csv")
df <- read.csv(file_path)

num_cols <- c("threshold", "predicted_correlation",
              "SynGWAS_hat_beta", "SynGWAS_sd_hat_beta", "SynGWAS_pvalues",
              "ObsGWAS_hat_beta", "ObsGWAS_sd_hat_beta", "ObsGWAS_pvalues")

df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(as.character(x)))

# ==========================================
# 2. 计算染色体全局偏移量
# ==========================================
chr_order <- c(as.character(1:22))

chr_offsets <- final_protein_map %>%
  filter(chr %in% chr_order) %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos)) %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  arrange(chr) %>%
  mutate(offset = cumsum(as.numeric(lag(chr_len, default = 0))))

# ==========================================
# 3. 核心作图数据对齐 (plot_df)
# ==========================================
df$chr <- as.character(df$chr)
chr_offsets$chr <- as.character(chr_offsets$chr)
final_protein_map$chr <- as.character(final_protein_map$chr)

## revision: the Syn-PALM primary associations at P < 5e-8 (as in the text), not the
## union of Syn-PALM and Obs-GWAS lead variants at the second threshold (1e-8)
df_sub <- as.data.frame(data.table::fread(
  "<PROJECT_DIR>/results/table_results_cf_noPRS_noPC/loci_SynG_5e-08.csv")[
  , .(chr = as.character(chr), POS38, rsid, protein_short_name = gsub("-", "_", protein),
      SynGWAS_pvalues = 10^(-SynSurrG_negative_log10_pvalues))])

plot_df <- df_sub %>%
  inner_join(chr_offsets, by = "chr") %>%
  mutate(qtl_pos_cum = as.numeric(POS38) + as.numeric(offset)) %>%
  left_join(final_protein_map %>% select(id, gene_chr = chr, gene_pos = pos), 
            by = c("protein_short_name" = "id")) %>%
  left_join(chr_offsets, by = c("gene_chr" = "chr"), suffix = c(".qtl", ".gene")) %>%
  mutate(gene_pos_cum = as.numeric(gene_pos) + as.numeric(offset.gene)) %>%
  mutate(is_cis = if_else(chr == gene_chr & abs(as.numeric(POS38) - as.numeric(gene_pos)) < 1e6, 
                          "Cis", "Trans"))

## revision: hotspots are labelled by the nearest protein-coding gene (TxDb hg38) of the most
## frequent lead variant in the window, not by a protein name (each protein
## contributes at most a few lead variants to a window, so the most frequent
## protein name was effectively the first in alphabetical order)
source("<PROJECT_DIR>/code/post_cf3/nearest_coding_gene.R")
.ng <- unique(plot_df[, c("rsid", "chr", "POS38")])
.ng <- nearest_coding_gene(.ng$rsid, .ng$chr, .ng$POS38)
rsid2gene <- setNames(.ng$nearest_gene, .ng$rsid)
top_gene <- function(r) { g <- rsid2gene[names(sort(table(r), decreasing = TRUE))[1]]; if (is.na(g)) names(sort(table(r), decreasing = TRUE))[1] else unname(g) }

# ==========================================
# 4. 生成中心主图 (p_main)
# ==========================================
## revision: chromosomes 18, 20 and 21 are not labelled (the labels of the short chromosomes ran together)
chr_lab3 <- ifelse(chr_offsets$chr %in% c("18", "20", "21"), "", chr_offsets$chr)
p_main <- ggplot(data = plot_df, 
                 aes(x = qtl_pos_cum, y = gene_pos_cum)) +
  geom_vline(xintercept = chr_offsets$offset, color = "grey90") +
  geom_hline(yintercept = chr_offsets$offset, color = "grey90") +
  geom_point(aes(color = is_cis, size = -log10(SynGWAS_pvalues)), alpha = 0.4) +
  scale_color_manual(values = c("Cis" = "#E15759", "Trans" = "#76B7B2"), 
                     guide = guide_legend(title = NULL, override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.01, 2)) +
  scale_x_continuous(breaks = chr_offsets$offset + chr_offsets$chr_len/2, labels = chr_lab3, expand = c(0,0)) +
  scale_y_continuous(breaks = chr_offsets$offset + chr_offsets$chr_len/2, labels = chr_lab3, expand = c(0,0)) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(x = "pQTL Position (Chromosome)", y = "Protein Position (Chromosome)", size = expression(-log[10] * P))

# ==========================================
# 5. 生成顶部热点图 (p_top)
# ==========================================
window_size <- 5e6 # 5MB 窗口

hotspot_labels <- plot_df %>%
  mutate(window = floor(qtl_pos_cum / window_size)) %>%
  group_by(window) %>%
  summarise(
    original_count = n(),
    top_protein = top_gene(rsid),
    pos_cum = mean(qtl_pos_cum, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(display_count = ifelse(original_count > 250, 250, original_count)) %>%
  mutate(label_text = ifelse(original_count > 250, 
                             paste0(top_protein, " (", original_count, ")"), 
                             top_protein)) %>%
  arrange(desc(original_count)) %>%
  slice_head(n = 11)

p_top <- ggplot(plot_df %>% 
                  mutate(window = floor(qtl_pos_cum / window_size)) %>%
                  group_by(window) %>% 
                  summarise(count = pmin(n(), 250), pos = mean(qtl_pos_cum), .groups = 'drop'), 
                aes(x = pos, y = count)) +
  geom_col(fill = "grey40", width = window_size) + 
  geom_text_repel(
    data = hotspot_labels,
    aes(x = pos_cum, y = display_count, label = label_text),
    size = 3.2, fontface = "italic", direction = "both", nudge_y = 110, box.padding = 0.25,
    segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf, force = 8, seed = 1
  ) +
  scale_y_continuous(limits = c(0, 560), breaks = c(0, 100, 200), expand = c(0, 0), guide = guide_axis(cap = "upper")) + 
  scale_x_continuous(expand = c(0,0), limits = c(0, max(chr_offsets$offset + chr_offsets$chr_len) + 1e7)) +
  theme_classic(base_size = 13) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), plot.margin = margin(b = 0, t = 20, r = 10, l = 10)) +
  labs(y = "No. of associated\nproteins")

# ==========================================
# 6. 生成右侧变异数量统计图 (p_right)
# ==========================================
protein_counts <- plot_df %>%
  group_by(protein_short_name) %>%
  summarise(
    count = n(),
    gene_pos_cum = unique(gene_pos_cum),
    .groups = 'drop'
  ) %>%
  mutate(display_count = pmin(count, 50)) # 封顶值，防止个别蛋白关联太多变异导致图表失衡

p_right <- ggplot(protein_counts, aes(y = gene_pos_cum, x = display_count)) +
  geom_col(orientation = "y", fill = "grey40", width = window_size) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) + 
  # 重要：这里的 Y 轴 limits 必须与中心主图保持绝对一致
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(chr_offsets$offset + chr_offsets$chr_len) + 1e7)) +
  theme_classic(base_size = 13) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(l = 0, r = 20, t = 10, b = 10)) +
  labs(x = "No. of \nsignificant variants")

# ==========================================
# 7. 拼接与导出 (Patchwork)
# ==========================================
combined_plot <- p_top + plot_spacer() + p_main + p_right + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1), 
    heights = c(1.5, 4),
    guides = 'collect'
  ) & theme(legend.position = "bottom")

ggsave(
  plot = combined_plot,
  file = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/Summary_of_pQTLs_ThreePanel.pdf",
  width = 8.5,    # revision: smaller canvas, larger text (panel A of Fig 5 is about 12.5 cm wide)
  height = 7.6, 
  units = "in",
  dpi = 480
)

# ==========================================
# 5. 生成顶部热点图 (p_top) — 优化版
# ==========================================

window_size <- 5e6

# --- 1. 汇总窗口数据 ---
window_df <- plot_df %>%
  mutate(window = floor(qtl_pos_cum / window_size)) %>%
  group_by(window) %>%
  summarise(
    original_count  = n(),
    top_protein     = top_gene(rsid),
    pos_cum         = mean(qtl_pos_cum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    display_count = pmin(original_count, 250),
    is_hotspot    = original_count > 250
  )

# --- 2. 筛选标签（count > 250 的前 11 个） ---
hotspot_labels <- window_df %>%
  filter(original_count > 250) %>%
  arrange(desc(original_count)) %>%
  slice_head(n = 11) %>%
  mutate(label_text = paste0(top_protein, " (", original_count, ")"))

# --- 3. 配色方案 ---
col_hotspot <- "#94a3b8"
col_normal  <- "#94a3b8"
col_label   <- "grey10"
col_segment <- "grey10"

# --- 4. 绘图 ---
p_top <- ggplot(window_df, aes(x = pos_cum, y = display_count)) +
  
  geom_col(
    aes(fill = is_hotspot),
    width = window_size * 0.92,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("FALSE" = col_normal, "TRUE" = col_hotspot)) +
  
  geom_col(
    data = filter(window_df, is_hotspot),
    aes(x = pos_cum, y = display_count),
    fill = col_hotspot,
    width = window_size * 0.92,
    color = "#94a3b8",
    linewidth = 0.25
  ) +
  
  geom_text_repel(
    data             = hotspot_labels,
    aes(x = pos_cum, y = display_count, label = label_text),
    size             = 2.8,
    fontface         = "italic",
    color            = col_label,
    direction        = "y",
    nudge_y          = 28,
    segment.color    = col_segment,
    segment.size     = 0.3,
    segment.linetype = "dashed",
    min.segment.length = 0,
    max.overlaps     = Inf,
    force            = 6,
    box.padding      = 0.2,
    point.padding    = 0.1
  ) +
  
  scale_y_continuous(
    limits = c(0, 300),
    breaks = c(0, 100, 200, 300),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = chr_offsets %>%
      mutate(mid = offset + chr_len / 2) %>%
      pull(mid),
    labels = chr_offsets$chr,
    expand = c(0, 0),
    limits = c(0, max(chr_offsets$offset + chr_offsets$chr_len) + 1e7)
  ) +
  
  geom_vline(
    xintercept = chr_offsets$offset[-1],
    color      = "grey80",
    linewidth  = 0.25,
    linetype   = "dashed"
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    # x 轴 — 加粗
    axis.title.x = element_text(size = 9,   color = "grey10", face = "bold",
                                margin = margin(t = 6)),
    axis.text.x  = element_text(size = 7.5, color = "grey10", face = "bold"),
    axis.ticks.x   = element_line(color = "grey40", linewidth = 0.6),
    axis.ticks.length.x = unit(4, "pt"),
    axis.line.x    = element_line(color = "grey40", linewidth = 0.6),
    
    # y 轴 — 加粗
    axis.title.y = element_text(size = 9,   color = "grey10", face = "bold",
                                lineheight = 1.3),
    axis.text.y  = element_text(size = 8,   color = "grey10", face = "bold"),
    axis.ticks.y   = element_line(color = "grey40", linewidth = 0.6),
    axis.ticks.length.y = unit(4, "pt"),
    axis.line.y    = element_line(color = "grey40", linewidth = 0.6),
    
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 24, r = 10, b = 4, l = 10)
  ) +
  
  labs(
    y = "No. of associated proteins",
    x = "pQTL Position (Chromosome)"
  )

p_top

ggsave(
  filename = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/p_top_hotspot.pdf",
  plot     = p_top,
  width    = 8,
  height   = 4,
  units    = "in",
  device   = cairo_pdf
)

##########

library(ggplot2)
library(dplyr)

window_size <- 5e6

# --- 1. 汇总窗口数据 ---
window_df <- plot_df %>%
  mutate(window = floor(qtl_pos_cum / window_size)) %>%
  group_by(window) %>%
  summarise(
    original_count = n(),
    top_protein    = top_gene(rsid),
    pos_cum        = mean(qtl_pos_cum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    display_count = pmin(original_count, 250),
    is_hotspot    = original_count > 250
  )

# --- 2. 染色体信息 ---
chr_info <- plot_df %>%
  group_by(chr) %>%
  summarise(
    start_cum = min(qtl_pos_cum),
    end_cum   = max(qtl_pos_cum),
    mid_cum   = (min(qtl_pos_cum) + max(qtl_pos_cum)) / 2,
    .groups   = "drop"
  ) %>%
  mutate(chr_num = as.integer(as.character(chr))) %>%
  arrange(chr_num)

chr_bounds <- chr_info %>% filter(chr_num > 1)
genome_max <- max(plot_df$qtl_pos_cum)

# --- 3. hotspot 标签 + 计算旋转角度 ---
hotspot_labels <- window_df %>%
  filter(is_hotspot) %>%
  arrange(desc(original_count)) %>%
  slice_head(n = 11) %>%
  ## revision: hotspots closer than 2.5% of the genome would print on top of each
  ## other; they share one label (most frequent first), placed at their mean position
  arrange(pos_cum) %>%
  mutate(grp = cumsum(c(1, diff(pos_cum) > 0.025 * max(plot_df$qtl_pos_cum)))) %>%
  group_by(grp) %>% arrange(desc(original_count), .by_group = TRUE) %>%
  summarise(label_one = paste0(top_protein, " (", original_count, ")", collapse = "\n"),
            pos_cum = mean(pos_cum), original_count = max(original_count), .groups = "drop") %>%
  mutate(
    label_text = ifelse(grepl("\n", label_one), label_one, sub(" \\(", "\n(", label_one)),
    # 角度计算：coord_polar start=0 从顶部顺时针
    angle_deg  = 90 - (pos_cum / genome_max) * 360,
    # 左半圈（pos_cum > genome_max/2）翻转 180° 保持可读
    angle_text = ifelse(pos_cum > genome_max / 2,
                        angle_deg + 180,
                        angle_deg),
    # 左侧右对齐，右侧左对齐
    hjust_val  = ifelse(pos_cum > genome_max / 2, 1, 0)
  )

# --- 4. y 轴刻度 ---
y_tick_x <- min(plot_df$qtl_pos_cum)
y_ticks  <- c(0, 100, 200)

DATA_TOP <- 255

# --- 5. 绘图 ---
p_circular <- ggplot() +
  
  # 染色体分隔线（淡）
  geom_vline(
    data = chr_bounds,
    aes(xintercept = start_cum),
    color = "grey90", linewidth = 0.15
  ) +
  
  # 数据柱子
  geom_col(
    data  = window_df,
    aes(x = pos_cum, y = display_count, fill = is_hotspot),
    width = window_size * 0.90,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("FALSE" = "#94a3b8", "TRUE" = "#94a3b8")) +
  
  # y 轴刻度标签
  annotate("text",
           x = y_tick_x, y = y_ticks,
           label = as.character(y_ticks),
           size = 2.3, color = "grey40", hjust = 1.3
  ) +
  annotate("segment",
           x    = y_tick_x - genome_max * 0.003, xend = y_tick_x,
           y    = y_ticks, yend = y_ticks,
           color = "grey55", linewidth = 0.3
  ) +
  
  # ★ 蛋白质标签：辐射状向外发散
  geom_text(
    data = hotspot_labels,
    aes(x     = pos_cum,
        y     = DATA_TOP + 10,
        label = label_text,
        angle = angle_text,
        hjust = hjust_val),
    size       = 3.8,
    fontface   = "bold",
    color      = "grey10",
    lineheight = 0.85,
    vjust      = 0.5
  ) +
  
  # ★ 染色体号：圆心内圈，浅蓝色
  geom_text(
    data = chr_info,
    aes(x = mid_cum, y = -10, label = chr_num),
    size     = 2,
    color    = "#7eb8d4",
    fontface = "bold"
  ) +
  
  # 坐标轴
  scale_x_continuous(
    limits = c(0, genome_max),
    breaks = chr_info$mid_cum,
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-150, DATA_TOP + 10),
    expand = c(0, 0)
  ) +
  
  coord_polar(theta = "x", start = 0, clip = "off") +
  
  theme_minimal(base_size = 10) +
  theme(
    axis.line          = element_blank(),
    axis.ticks         = element_blank(),
    axis.text          = element_blank(),
    axis.title         = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color     = "grey92",
                                      linewidth = 0.25,
                                      linetype  = "dashed"),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(30, 30, 30, 30)
  ) +
  theme(
    plot.caption = element_text(size = 7, color = "grey50",
                                hjust = 0.5, margin = margin(t = 4))
  )

p_circular

ggsave(
  filename = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/p_top_hotspot_cr.pdf",
  plot     = p_circular,
  width    = 8,
  height   = 8,
  units    = "in",
  device   = cairo_pdf
)


# ==========================================
# p_right — 蛋白关联变异数柱状图（优化版，方向同 p_top）
# ==========================================

col_hotspot <- "#94a3b8"
col_normal  <- "#94a3b8"

protein_counts <- plot_df %>%
  group_by(protein_short_name) %>%
  summarise(
    count        = n(),
    gene_pos_cum = unique(gene_pos_cum),
    .groups      = "drop"
  ) %>%
  mutate(
    display_count = pmin(count, 50),
    is_hotspot    = count > 20
  )

p_right <- ggplot(protein_counts, aes(x = gene_pos_cum, y = display_count)) +
  
  geom_col(
    aes(fill = is_hotspot),
    width       = window_size * 0.92,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("FALSE" = col_normal, "TRUE" = col_hotspot)) +
  
  geom_vline(
    xintercept = chr_offsets$offset[-1],
    color      = "grey80",
    linewidth  = 0.25,
    linetype   = "dashed"
  ) +
  
  scale_x_continuous(
    breaks = chr_offsets %>%
      mutate(mid = offset + chr_len / 2) %>%
      pull(mid),
    labels = chr_offsets$chr,
    expand = c(0, 0),
    limits = c(0, max(chr_offsets$offset + chr_offsets$chr_len) + 1e7)
  ) +
  
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 60),
    breaks = c(0, 20, 40, 60)
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    # x 轴 — 加粗
    axis.title.x = element_text(size = 9,   color = "grey10", face = "bold",
                                margin = margin(t = 6)),
    axis.text.x  = element_text(size = 7.5, color = "grey10", face = "bold"),
    axis.ticks.x   = element_line(color = "grey40", linewidth = 0.6),
    axis.ticks.length.x = unit(4, "pt"),
    axis.line.x    = element_line(color = "grey40", linewidth = 0.6),
    
    # y 轴 — 加粗
    axis.title.y = element_text(size = 9,   color = "grey10", face = "bold",
                                lineheight = 1.3),
    axis.text.y  = element_text(size = 8,   color = "grey10", face = "bold"),
    axis.ticks.y   = element_line(color = "grey40", linewidth = 0.6),
    axis.ticks.length.y = unit(4, "pt"),
    axis.line.y    = element_line(color = "grey40", linewidth = 0.6),
    
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 4, r = 10, b = 10, l = 10)
  ) +
  
  labs(
    x = "Protein Position (Chromosome)",
    y = "No. of significant variants"
  )

p_right

ggsave(
  filename = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/p_right_hotspot.pdf",
  plot     = p_right,
  width    = 8,
  height   = 4,
  units    = "in",
  device   = cairo_pdf
)

library(ggplot2)
library(dplyr)

# --- 1. 蛋白质数据 ---
protein_counts <- plot_df %>%
  group_by(protein_short_name) %>%
  summarise(
    count        = n(),
    gene_pos_cum = mean(gene_pos_cum, na.rm = TRUE),
    top_rsid     = rsid[which.min(SynGWAS_pvalues)],
    .groups      = "drop"
  ) %>%
  mutate(
    display_count = pmin(count, 50),
    is_hotspot    = display_count > 45
  )

cat(sprintf("蛋白质总数: %d，hotspot: %d\n",
            nrow(protein_counts), sum(protein_counts$is_hotspot)))

# --- 2. 染色体边界 ---
x_max <- max(chr_offsets$offset + chr_offsets$chr_len) + 1e7
genome_max <- x_max

chr_bounds <- data.frame(xintercept = chr_offsets$offset[-1])
chr_mids   <- chr_offsets %>% mutate(mid = offset + chr_len / 2)

# --- 3. labels: the N_LAB proteins with the most primary associations (revision 2026-10-19;
##     the earlier rule, more than 45 associations, labelled only 2 proteins). Proteins whose
##     genes lie within 1.2% of the genome of each other share one label, so labels do not overlap.
N_LAB <- 20
top_prot <- protein_counts %>% filter(!is.na(gene_pos_cum)) %>%
  slice_max(count, n = N_LAB, with_ties = FALSE)
protein_counts <- protein_counts %>% mutate(is_hotspot = protein_short_name %in% top_prot$protein_short_name)
hotspot_proteins <- top_prot %>%
  arrange(gene_pos_cum) %>%
  mutate(grp = cumsum(c(TRUE, diff(gene_pos_cum) > genome_max * 0.012))) %>%
  group_by(grp) %>% arrange(desc(count), .by_group = TRUE) %>%
  summarise(label_text   = paste0(toupper(protein_short_name), " (", count, ")", collapse = "\n"),
            gene_pos_cum = mean(gene_pos_cum), .groups = "drop") %>%
  mutate(
    angle_deg  = 90 - (gene_pos_cum / genome_max) * 360,
    angle_text = ifelse(gene_pos_cum > genome_max / 2, angle_deg + 180, angle_deg),
    hjust_val  = ifelse(gene_pos_cum > genome_max / 2, 1, 0)
  )

cat(sprintf("标注蛋白质数: %d\n", nrow(hotspot_proteins)))

# --- 4. 参数 ---
y_tick_x    <- x_max * 0.995
y_ticks     <- c(0, 20, 40)
DATA_TOP    <- 50
HOLE_SIZE   <- -30
CHR_LABEL_Y <- -3

# --- 5. 绘图 ---
p_right_circular <- ggplot() +
  
  # 染色体分隔线（细实线）
  geom_segment(
    data = chr_bounds,
    aes(x = xintercept, xend = xintercept,
        y = HOLE_SIZE, yend = DATA_TOP),
    color     = "grey85",
    linewidth = 0.20,
    linetype  = "solid"
  ) +
  
  # 数据柱子
  geom_col(
    data  = protein_counts,
    aes(x = gene_pos_cum, y = display_count, fill = is_hotspot),
    width = genome_max / 1500,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("FALSE" = "#b6c2d1", "TRUE" = "#334155")) +   # labelled proteins darker
  
  # ★ y 轴数字（无横杠，居中）
  annotate("text",
           x     = y_tick_x,
           y     = y_ticks,
           label = as.character(y_ticks),
           size  = 2.3, color = "grey50", hjust = 0.5
  ) +
  
  # 蛋白质标签：辐射状向外
  geom_text(
    data = hotspot_proteins,
    aes(x     = gene_pos_cum,
        y     = DATA_TOP + 5,
        label = label_text,
        angle = angle_text,
        hjust = hjust_val),
    size       = 3.6,
    fontface   = "bold",
    color      = "grey10",
    lineheight = 0.85,
    vjust      = 0.5
  ) +
  
  # 染色体号：内圈，浅蓝色
  geom_text(
    data = chr_mids,
    aes(x = mid, y = CHR_LABEL_Y, label = chr),
    size     = 2.0,
    color    = "#7eb8d4",
    fontface = "bold"
  ) +
  
  # 坐标轴
  scale_x_continuous(
    limits = c(0, genome_max),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(HOLE_SIZE, DATA_TOP + 10),
    breaks = y_ticks,
    expand = c(0, 0)
  ) +
  
  coord_polar(theta = "x", start = 0, clip = "off") +
  
  theme_minimal(base_size = 10) +
  theme(
    axis.line          = element_blank(),
    axis.ticks         = element_blank(),
    axis.text          = element_blank(),
    axis.title         = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    # ★ 同心虚线圆作为 y 轴参考
    panel.grid.major.y = element_line(color     = "grey88",
                                      linewidth = 0.25,
                                      linetype  = "dashed"),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(30, 30, 30, 30)
  ) +
  theme(
    plot.caption = element_text(size = 7, color = "grey50",
                                hjust = 0.5, margin = margin(t = 4))
  )

p_right_circular

ggsave(
  filename = "<PROJECT_DIR>/results/figs_cf_noPRS_noPC/p_right_hotspot_cr.pdf",
  plot     = p_right_circular,
  width    = 8,
  height   = 8,
  units    = "in",
  device   = cairo_pdf
)
