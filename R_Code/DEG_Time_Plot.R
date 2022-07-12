common_des <- rownames(read.csv('./../gPD_HC/outputData/common_de.csv' , row.names = 1))

deg_all_time <- c()
for (TS in c('BL' , 'V02' , 'V04' , 'V06' , 'V08')) {
  load(paste0('./../gPD_HC/preprocessedData/',TS,'/preprocessed2.RData'))
  tmp <- datExpr[common_des , ] %>% t %>% data.frame %>% 
    rownames_to_column('ID') %>% left_join(datMeta %>% rownames_to_column('ID') , by = 'ID') %>% 
    select(ID, common_des , meta_df.Clinical_Event , PD_status)
  deg_all_time <- rbind(deg_all_time , tmp)
}

deg_all_time <- datExpr[common_des[1] , ] %>% data.frame %>% setNames(. , common_des[1])


ggplot(deg_all_time , aes(x = meta_df.Clinical_Event , y = ENSG00000221869.4  , fill = PD_status)) + geom_boxplot() 
       