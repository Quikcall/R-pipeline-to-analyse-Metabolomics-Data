#1H NMR Metabolomics Identifies Underlying Inflammatory Pathology in Osteoarthritis and Rheumatoid Arthritis Synovial Joints.



#1 - Curation and pre-processing of data
#Data Reading into specmine
get_metabolights_study_files_assay("MTBLS564", 1, "C:/project_extraction/MTBLS564") #Ficheiros de dados
get_metabolights_study_metadata_assay("MTBLS564", 1, "C:/project_extraction/MTBLS564") #Ficheiro de metadados
get_metabolights_study_samples_files("MTBLS564", 1, "C:/project_extraction/MTBLS564") #Ficheiro samples_files, que faz a correspondencia



library(caret)
library(specmine)
read_csvs_folder("C:/project_extraction/MTBLS564", header=TRUE, sep=",", dec=".")


                                                                                      #entre os nomes das amostras e os nomes dos ficheiros
MTBLS564_data=read_dataset_csv(filename.data="C:/project_extraction/MTBLS564/hSF_MTBLS_normalised_2.csv",
filename.meta="C:/project_extraction/MTBLS564/metadata1.csv", type="concentrations", description =
                                  "Metabolite profiling from MetaboLights study MTBLS564", label.x = "Metabolites", label.values =
                                  "intensity")

#data visualization
sum_dataset(MTBLS564_data)
MTBLS564_data$data[1:6,1:6]
MTBLS564_data$metadata[1:6,]

#data subsetting
MTBLS564_OA_RA=subset_samples_by_metadata_values(MTBLS564_data,"Diagnosis",c("OA", "RA"))

#2 - DESCRIPTIVE STATISTICS 
#t test

ttest_MTBLS564_OA_RA=tTests_dataset(MTBLS564_OA_RA, "Diagnosis")
#t test result for significant p values
ttest_MTBLS564_OA_RA[ttest_MTBLS564_OA_RA$fdr<0.05,]
significant_MTBLS564 = ttest_MTBLS564_OA_RA[ttest_MTBLS564_OA_RA$fdr<0.05,]



#select metabolites that are different between diagnosis
#Metabolites that will remain:
metabolites = rownames(significant_MTBLS564)
#Taking the significant metabolites from the diagnosis aside:
MTBLS564_nrm_metabolites =  subset_x_values(MTBLS564_data, metabolites)
#Get a dataset with only the samples from OA and RA diagnosis:
MTBLS564_OA_RA_new = subset_samples_by_metadata_values(MTBLS564_nrm_metabolites, "Diagnosis", c("OA", "RA"))

#Remove unknown samples
significant_MTBLS564_known = subset_x_values(MTBLS564_OA_RA_new, grep(".*[U,u]nknown.*", rownames(MTBLS564_OA_RA_new$data), value=T, invert=T))


#ttest for significant and known metabolites
ttest_significant_MTBLS564_known = tTests_dataset(significant_MTBLS564_known, "Diagnosis")
#plot of ttest (pvalues) with threshold of 0.05
plot_ttests(significant_MTBLS564_known, ttest_significant_MTBLS564_known, tt.threshold = 0.05)

#Anova 1 way 
anova_res=aov_all_vars(significant_MTBLS564_known, "Diagnosis", doTukey = T, write.file = F, file.out = "anova-res.csv")
plot_anova(significant_MTBLS564_known, anova_res, anova.threshold = 0.05, reverse.x = F)
#confirms the significative results from the t test

#Kruskal-Wallis test
kruskal_res=kruskalTest_dataset(significant_MTBLS564_known, "Diagnosis", threshold = 0.01, write.file = F, file.out = "kruskal.csv") 
plot_kruskaltest(significant_MTBLS564_known, kruskal_res, kr.threshold = 0.01)


##DIFFERENTIAL EXPRESSION
#boxplot of significant metabolites
boxplot_variables(significant_MTBLS564_known, horizontal=F, nhar.label=100, cex.axis=0.5)
boxplot_vars_factor(significant_MTBLS564_known, "Diagnosis", variables = row.names(significant_MTBLS564_known$data), cex.axis=1, cex.main=1, vec.par=c(3,5))
#meter nomes únicos caso seja de interesse
#make.names(metabolites$metabolite_identification, unique=TRUE)


#fold change
fc_MTBLS564_OA_known = fold_change(significant_MTBLS564_known, "Diagnosis", "OA")
fc_MTBLS564_RA_known = fold_change(significant_MTBLS564_known, "Diagnosis", "RA")
#fold change plot
par(mfrow = c(1, 2))
plot_fold_change(significant_MTBLS564_known, na.omit(fc_MTBLS564_OA_known), 2)
plot_fold_change(significant_MTBLS564_known, na.omit(fc_MTBLS564_RA_known), 2)
par(mfrow=c(1, 1))


#volcano plot
volcano_plot_fc_tt(significant_MTBLS564_known, na.omit(fc_MTBLS564_OA_known), ttest_significant_MTBLS564_known[rownames(fc_MTBLS564_OA_known),], fc.threshold=2, tt.threshold = 0.05)
volcano_plot_fc_tt(significant_MTBLS564_known, na.omit(fc_MTBLS564_RA_known), ttest_significant_MTBLS564_known[rownames(fc_MTBLS564_RA_known),], fc.threshold=2, tt.threshold = 0.05)


#Heatmap
colors = rev(colorRampPalette(RColorBrewer::brewer.pal(10,"RdBu"))(526))
stats::heatmap(MTBLS564_OA_RA$data, scale="none", Colv=NA, margins=c(0,8), cexRow=0.6, cexCol=.4, keep.dendro=FALSE, col=colors, labCol=rep("", 34))
stats::heatmap(significant_MTBLS564_known$data, scale="none", Colv=significant_MTBLS564_known$metadata$Diagnosis, margins=c(0,8), cexRow=0.6, cexCol=.4, keep.dendro=TRUE, col=colors, rownames(significant_MTBLS564_known$metadata$Diagnosis))


####################################################################################################################################
##Enrichment Analysis & Pathway Analysis


#BiocManager::install("FELLA", version = "3.8")
library(FELLA)
signif_mets = row.names(significant_MTBLS564_known$data)
og_data = read.table('m_mtbls564_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv', header=TRUE, sep='\t')
all_met_ids = og_data[c(1,5,7)]
ids = c()
met_names = c()
all_met_ids$database_identifier = as.character(all_met_ids$database_identifier)
all_met_ids$metabolite_identification = as.character(all_met_ids$metabolite_identification)

for (i in 1:length(all_met_ids$pattern_ID)){
  if (is.element(all_met_ids$pattern_ID[i], signif_mets)==TRUE){ids[i]= all_met_ids$database_identifier[i];met_names[i]=all_met_ids$metabolite_identification[i]}}


met_ids = na.omit(ids)
met_names = unique(met_names) #excluindo as variantes, existem 32 metabolitos significativos, sendo q 2 deles são desconhecidos
met_names = na.omit(met_names)
met_names = unique(met_names)
#KEGG ids converted from CHEBI using MetaboAnalyst
met_db_ids = read.table('name_map.csv',header=TRUE,sep = ',')
met_kegg_id = data.frame(met_db_ids$KEGG)
met_kegg_id= met_kegg_id[!apply(met_kegg_id == "", 1, all),] #remove empty cells



#build database
set.seed(1)
# Filter overview pathways
graph <- buildGraphFromKEGGREST(  organism = "hsa", filter.path = c("01100", "01200", "01210", "01212", "01230"))
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph, databaseDir = tmpdir,internalDir = FALSE,matrices = "diffusion", normality = "diffusion",niter = 50)
fella.data <- loadKEGGdata( databaseDir = tmpdir, internalDir = FALSE, loadMatrix = "diffusion" )
fella.data
cat(getInfo(fella.data))

synovial_fluid_analysis = defineCompounds(compounds = as.character(met_kegg_id), data = fella.data) 
getInput(synovial_fluid_analysis)


#execute enrichment
#Enriching using diffusion
synovial_fluid_analysis = runDiffusion(object = synovial_fluid_analysis, data = fella.data, approx = "normality")

nlimit <- 50
vertex.label.cex <- 0.5
#plot pathways
plot(synovial_fluid_analysis, method = "diffusion", data = fella.data, nlimit = nlimit,vertex.label.cex = vertex.label.cex)

table_all <- generateResultsTable( method = "diffusion", nlimit = nlimit, object = synovial_fluid_analysis, data = fella.data)






################################################################################################################################################################################### 



MTBLS564_bruker = read_Bruker_files("C:/project_extraction/MTBLS564/1/data", metadata_file =
                                      "C:/project_extraction/MTBLS564/metadata1.csv", samples.names = "C:/project_extraction/MTBLS564/samples_files.csv",
                                    zipped=T, description="1H NMR Metabolomics", label.x="ppm", label.values="intensity")

sum_dataset(MTBLS564_bruker)
count_missing_values(MTBLS564_bruker) #116
stats_by_sample(MTBLS564_bruker, samples = NULL)
plot_spectra(MTBLS564_data, "SF spectra", cex=1)


#Missing values
mtb546_nona =missingvalues_imputation(MTBLS564_bruker, method = "value", value = 5e-4)

#Spectra plot
plot_spectra(mtb546_nona, 'Diagnosis', main="Spectra plot", variable.bounds = c(0,6), lty = 1, cex=1.2)
#Peak Detection
MTBLS546_peaks = detect_nmr_peaks_from_dataset(mtb546_nona, baseline_tresh = 50000)
plot_peaks(MTBLS546_peaks, "Diagnosis")

#Boxplot of variables along metadata variables:
boxplot_vars_factor(MTBLS546_peaks, "Diagnosis", variables = row.names(MTBLS546_peaks$data), cex.axis=0.5, cex.main=1, vec.par=c(3,3))
#Boxplot of variables distribution:
boxplot_vars_factor(MTBLS546_peaks, "Diagnosis", variables = c("3.05", "4.27"))



################################################################################################################################################################################### 

#Redução da Dimensionalidade
#	dado	um	conjunto	alargado	de	variáveis,	descobrir	um conjunto	mais	pequeno	de	variáveis	não	correlacionadas	entre	si	
#que	explicam	a	maior	parte	da	variabilidade	dos	dados usando a técnica de Análise de Componentes Principais (PCA) usando
#os dados significativos

pca_mtbls564 = pca_analysis_dataset(significant_MTBLS564_known, scale =F, center = F, write.file = F, file.out = "pca") #usa-se dataset normalizado

#importance table
library(caret)
library(ipred)
res_imp_pca_concentrations=pca_importance(pca_mtbls564, pcs=1:6) #coluna ordernadas por ordem decrescente dos valores de desvio padrão
DT::datatable(res_imp_pca_concentrations, options=list(scrollX=T))
#Plot showing results on all principal components
pca_screeplot(pca_mtbls564)


RColorBrewer::display.brewer.all()
#scores Plot
pca_scoresplot2D(significant_MTBLS564_known, pca_mtbls564, "Diagnosis", pcas = c(1,2), ellipses = TRUE, pallette = "Paired")
pca_scoresplot2D(significant_MTBLS564_known, pca_mtbls564, "Diagnosis", bw=TRUE)
pca_scoresplot3D(significant_MTBLS564_known, pca_mtbls564, "Diagnosis", pcas = c(1,2,3))
pca_scoresplot3D_rgl(significant_MTBLS564_known, pca_mtbls564, "Diagnosis", labels = TRUE)

#Biplot
pca_biplot(significant_MTBLS564_known, pca_mtbls564, colors="Diagnosis", cex = 0.0001, legend.cex = 0.5, x.colors = 1, inset = c(0, 0), legend.place = "topright")
pca_biplot3D(significant_MTBLS564_known, pca_mtbls564, "Diagnosis", pcas = c(1,2,3))

#Pairs Plot
pca_pairs_plot(significant_MTBLS564_known, pca_mtbls564, "Diagnosis", pcas = c(1,2,3))

#Kmeans Plot
pca_kmeans_plot2D(significant_MTBLS564_known, pca_mtbls564, num.clusters = 2)
pca_kmeans_plot3D(significant_MTBLS564_known, pca_mtbls564, num.clusters = 2, pcas = c(1,2,3))
#Kmeans Pairs Plot
pca_pairs_kmeans_plot(significant_MTBLS564_known, pca_mtbls564, num.clusters = 2, pcas = c(1,2,3,4,5))

#Cotovelo plot
wss = c()
for (k in 1:24) wss = c(wss, sum(kmeans(significant_MTBLS564_known$data, centers=k)$withinss) )
plot(1:24, wss, type="b", xlab="Num Clusters", ylab="WSS")


#############################################################################################################################################

#Clustering

cluster_mtbls564 = clustering(significant_MTBLS564_known, method = "hc", distance = "euclidean", type = "samples", num.clusters = 2, clustMethod = "complete")
#dendogram
dendrogram_plot_col(significant_MTBLS564_known, cluster_mtbls564, "Diagnosis", title = "Dendogram plot for MTBLS564 metabolites normalized dataset")
#K-Means Clustering
#Hierarchical clustering of the samples, using the euclidean distance and the agglomeration method complete:
res_kclust_concentrations = clustering(significant_MTBLS564_known, method="kmeans", num.clusters=2)
kmeans_plot(significant_MTBLS564_known, res_kclust_concentrations)






#############################################################################################################################################

#data = as.matrix(is.isignificant_MTBLS564_known$data)
#significant_MTBLS564_known$data = t(significant_MTBLS564_known$data)

set.seed(110)
res_train_model=train_models_performance(significant_MTBLS564_known, c("pls", "rf", "nnet"), "Diagnosis", "cv", metric="Accuracy")
res_train_model$performance
res_predict_model_pls=predict_samples(res_train_model$final.models$pls, significant_MTBLS564_known$data)
res_predict_model_rf=predict_samples(res_train_model$final.models$rf, significant_MTBLS564_known$data)
res_predict_model_nnet=predict_samples(res_train_model$final.models$nnet, significant_MTBLS564_known$data)


set.seed(111)
#res_train_model1=feature_selection(significant_MTBLS564_known, "Diagnosis", method="rfe", functions = caret::ldaFuncs, validation = "repeatedcv", number = 3, subsets = 2^(1:6))
#res_train_model1

res_train_model2=feature_selection(significant_MTBLS564_known, "Diagnosis", method="rfe", functions = caret::treebagFuncs, validation = "repeatedcv", number = 3, subsets = 2^(1:6))
res_train_model2 #resultados de fitting mais similares aos resultados de importância das variáveis do nosso modelo ao nosso modelo. O método é ideal para Regression e Classification models

#res_train_model3=feature_selection(significant_MTBLS564_known, "Diagnosis", method="rfe", functions = caret::rfFuncs, validation = "repeatedcv", number = 3, subsets = 2^(1:6))
#res_train_model3

var_imp1 = varImp(res_train_model$final.models$pls)
Overall_Importance = var_imp1[order(var_imp1$Overall, decreasing = T),c(0,1)]
var_imp_ord= as.data.frame(Overall_Importance, row.names = rownames(var_imp1))
head(var_imp_ord,10)

var_imp2 = varImp(res_train_model$final.models$rf)
Overall_Importance = var_imp2[order(var_imp2$Overall, decreasing = T),c(0,1)]
var_imp2_ord= as.data.frame(Overall_Importance, row.names = rownames(var_imp2))
head(var_imp2_ord,10)

var_imp3 = varImp(res_train_model$final.models$nnet)
Overall_Importance = var_imp3[order(var_imp3$Overall, decreasing = T),c(0,1)]
var_imp3_ord= as.data.frame(Overall_Importance, row.names = rownames(var_imp1))
head(var_imp3_ord,10)


