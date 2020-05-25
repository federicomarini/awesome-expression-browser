# awesome-expression-browser

[![Build Status](https://travis-ci.org/federicomarini/awesome-expression-browser.svg?branch=master)](https://travis-ci.org/federicomarini/awesome-expression-browser)

A curated list of software and resources for exploring and visualizing (browsing) expression data, but not only limited to that.
Credits for the backbone of the structure go to Sean Davis and his [awesome-single-cell](https://github.com/seandavi/awesome-single-cell) repository.

[Contributions welcome](https://github.com/federicomarini/awesome-expression-browser/blob/master/contributing.md)!

## Software list

- [AMP RA](https://immunogenomics.io/ampra), publicly available RNA-seq and CyTOF for human synovial tissue from patients with rheumatoid arthritis (RA) or osteoarthritis (OA), visualized with Shiny. [bioRxiv](https://doi.org/10.1101/351130)
- [AMP SLE](https://immunogenomics.io/ampsle), publicly available RNA-seq for human kidney biopsies from patients with systemic lupus erythematosus (SLE) patients, visualized with Shiny. [bioRxiv](https://doi.org/10.1101/363051)
- [adultHSPC10X](https://gottgens-lab.stemcells.cam.ac.uk/adultHSPC10X/) as a companion to https://doi.org/10.1182/blood-2017-12-821413
- [Allen Brain Atlases](http://portal.brain-map.org), Allen Brain Atlases and Data (from the Allen Institute) - for example, referred to in this recent publication on whole brain spatial transcriptomics (https://www.biorxiv.org/content/10.1101/784181v1)
- [ALS Spatiotemporal gene expression Atlas](https://als-st.nygenome.org), companion to the MS http://science.sciencemag.org/content/364/6435/89
- [ASAP](https://asap.epfl.ch), Automated Single-cell Analysis Pipeline (https://doi.org/10.1093/bioinformatics/btx337)
- [BioTuring Browser](http://bioturing.com/product/bbrowser), standalone browser from BioTuring
- [Blood RNAexpress Atlas](https://blueprint.haem.cam.ac.uk/bloodatlas/), composed of two separate browsers, for [mRNA](https://blueprint.haem.cam.ac.uk/mRNA/) and [miRNA](https://blueprint.haem.cam.ac.uk/micro/), with the same functionality to visualize expression values across samples from the BLUEPRINT project. Publication here: https://www.biorxiv.org/content/10.1101/764613v1?rss=1
- [Brain Cell RNA-seq Browser](http://jiaqianwulab.org/braincell/RNASeq.html) from the Jiaqian Wu Lab, related to publication in https://www.jneurosci.org/content/34/36/11929.short
- [Brain Immune Atlas](http://www.brainimmuneatlas.org/), a resource and visualization tool for assessing various single-cell RNA sequencing datasets that together capture the diversity of the brain immune compartment. Publication related: https://www.nature.com/articles/s41593-019-0393-4; source code: https://github.com/saeyslab/brainimmuneatlas/
- [Brain RNA-Seq](http://www.brainrnaseq.org), from the Barres Lab. Allows to explore different brain-related RNA-seq datasets
- [BRAIN-SAT](http://brainsat.eu), for interactive analysis of glia transcriptome studies - publication available here: https://www.biorxiv.org/content/10.1101/750000v1, successor to GOAD
- [cBioPortal](https://www.cbioportal.org), cBioPortal for Cancer Genomics
- [CellView](https://mbolisetty.shinyapps.io/CellView/), a Shiny app to visualize and explore single cell datasets. Source code available here: https://github.com/mohanbolisetty/CellView, described in the manuscript https://www.biorxiv.org/content/10.1101/123810v1. Currently not so clear what the expected format should be?
- [cellxgene](https://github.com/chanzuckerberg/cellxgene), deployed e.g. for the Tabula Muris (https://tabula-muris-senis.ds.czbiohub.org) - a simpler browser is available at https://tabula-muris.ds.czbiohub.org. Another example of this in action is http://mouse.retina.gofflab.org, for the Developing Mouse Retina https://doi.org/10.1016/j.neuron.2019.04.010
- [cerebro](https://github.com/romanhaa/Cerebro), Cell Report Browser, interactive application for exploration of scRNA-seq data. [bioRxiv](https://www.biorxiv.org/content/10.1101/631705v1) -> now published here: https://doi.org/10.1093/bioinformatics/btz877
- [ChickHeartAtlas](http://www.biomedatlas.org/ChickHeartAtlas/HH12/chick_heart_gene_expression.html), a 3D molecular atlas of the chick embryonic heart in support to the data of https://www.biorxiv.org/content/10.1101/609032v1
- [CoDEx](Cortical Development Expression Viewer), available at http://geschwindlab.dgsom.ucla.edu/pages/codexviewer (redirects to http://solo.bmap.ucla.edu/shiny/webapp/). Publication: https://doi.org/10.1016/j.neuron.2019.06.011
- [DropViz](http://dropviz.org), for exploring the Mouse Brain (Macosko & McCarroll Labs). Source code is at https://github.com/broadinstitute/dropviz, publication related to this is https://www.sciencedirect.com/science/article/pii/S0092867418309553
- [DRSC RNA-seq Explorer](https://www.flyrnai.org/tools/rna_seq/web/showProject/1/plot_coord=0/sample_id=all) Lab website: https://www.flyrnai.org/scRNA/ https://www.biorxiv.org/content/10.1101/410423v1 + 
- [DVEX](https://shiny.mdc-berlin.de/DVEX/) accompanying http://science.sciencemag.org/content/358/6360/194 (www.dvex.org)
- [epiviz](http://epiviz.cbcb.umd.edu), also accessible via epivizr (http://bioconductor.org/packages/release/bioc/html/epivizr.html)
- [Evo-devo mammalian organs](https://apps.kaessmannlab.org/evodevoapp/), a Shiny app enabling the exploration of individual gene expression profiles across organs, developmental stages and species in mammalians. Linked to the data presented in https://www.nature.com/articles/s41586-019-1338-5.
- [Expression Atlas](https://www.ebi.ac.uk/gxa/home) for bulk RNA-seq datasets
- [eyeIntegration](https://eyeIntegration.nei.nih.gov), all publicly available human eye RNA-seq datasets reprocessed and made available with a reactive Shiny app
- [FASTGenomics](https://fastgenomics.org), an online platform to share single-cell RNA sequencing data and analyses using reproducible workflows. Users can upload their own data and generate reproducible workflows ([Scholz et al. 2018](https://doi.org/10.1101/272476)). GitHub: [@fastgenomics](https://github.com/fastgenomics). A free demo instance is reachable at https://prod.fastgenomics.org/ (allows anonymous login)
- [Fate Bias Inference in Lymphoid progenitors](http://hematopoietic-progenitors.ie-freiburg.mpg.de), from the lab of Dominic Grün - single-cell RNA-seq data of murine hematopoietic progenitors
- [FireBrowse](http://firebrowse.org), developed by the Broad Institute on top of Firehose
- [Gene and protein expression in adult haematopoiesis](http://blood.stemcells.cam.ac.uk/single_cell_atlas.html), visualization of gene expression in 1,656 single HSPCs, accompanying https://doi.org/10.1182/blood-2016-05-716480
- [GenePattern](http://genepattern-notebook.org), Gene Pattern Notebook Environment (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5572818/)
- [GEPIA](http://gepia.cancer-pku.cn/) Gene Expression Profiling Interactive Analysis, a web server for cancer and normal gene expression profiling and interactive analyses. (https://academic.oup.com/nar/article/45/W1/W98/3605636)
- [Giotto](http://spatial.rc.fas.harvard.edu), a complete pipeline for integrative analysis and visualization of single-cell spatial transcriptomic data. Publication about it is https://doi.org/10.1101/701680 - Giotto Viewer is http://spatial.rc.fas.harvard.edu/giotto-viewer/, while the pipeline is available at https://github.com/RubD/Giotto
- [GTEx Portal](https://gtexportal.org/home/) Portal of the Genotype-Tissue Expression (GTEx) project 
- [Haemosphere](https://www.haemosphere.org/), to access and analyse the Haemopedia datasets and some other key datasets of interest to the haematopoiesis community (latest reference: https://doi.org/10.1093/nar/gky1020)
- [HCA Bone Marrow Viewer](http://www.altanalyze.org/ICGS/HCA/Viewer.php), related to this publication https://doi.org/10.1016/j.exphem.2018.09.004
- [HTSVis](http://htsvis.dkfz.de/HTSvis/), for image screening. Code available at https://github.com/boutroslab/HTSvis, described in https://academic.oup.com/bioinformatics/article/33/18/2960/3827333
- [Human Cell Atlas Data Portals](https://developmentcellatlas.ncl.ac.uk/datasets/hca_liver/), for exploring gene expression and developmental trajectories, showcased in this publication https://doi.org/10.1038/s41586-019-1652-y (HCA Liver instances: https://developmentcellatlas.ncl.ac.uk/datasets/hca_liver/). Source code is available at https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data
- [Human Liver Cell Atlas](http://human-liver-cell-atlas.ie-freiburg.mpg.de), related to https://dx.doi.org/10.1038/s41586-019-1373-2 - *currently not working?*
- [iMAC](https://www.stemformatics.org/atlas), an interactive atlas to explore phenotypic differences between macrophage cells. Puiblication available at https://www.biorxiv.org/content/10.1101/719237v1
- [indeXplorer](http://steinmetzlab.embl.de/shiny/indexplorer/), a web-based app for the explorative analysis of the hematopoietic stem cell dataset, accompanying https://doi.org/10.1038/ncb3493, with source code available at https://git.embl.de/velten/indeXplorer
- [Innate T cell gradient](https://immunogenomics.io/itc), low-input and single-cell RNA-seq for human T cell subsets, visualized with Shiny. [Nature Communications](https://doi.org/10.1038/s41467-019-08604-4)
- [IPF Cell Atlas](http://ipfcellatlas.com), the Idiopathic Pulmonary Fibrosis Cell Atlas. Linked to the publications https://www.biorxiv.org/content/10.1101/753806v1 and https://www.biorxiv.org/content/10.1101/759902v2
- [iS-CellR](https://github.com/immcore/iS-CellR), for analysing and visualising scRNA-seq data. R Shiny app. Manuscript available at https://academic.oup.com/bioinformatics/article/34/24/4305/5048937
- [iSEE](https://bioconductor.org/packages/release/bioc/html/iSEE.html) - Interactive SummarizedExperiment Explorer (dev version at https://github.com/iSEE/iSEE), with custom panels (https://github.com/iSEE/iSEE_custom). Publication available at https://f1000research.com/articles/7-741/v1
- [islet regulome](http://www.isletregulome.com/isletregulome/), to explore pancreatic islet genomic, epigenomic, and transcriptomic data. Data from https://www.nature.com/articles/s41588-019-0524-6, browser described in https://www.frontiersin.org/articles/10.3389/fgene.2017.00013/full. Source code is available at https://github.com/mireia-bioinfo/plotRegulome (plotting functions as a package!) and https://github.com/mireia-bioinfo/isletregulome_shiny (app)
- [iSTARTRAC](http://crctcell.cancer-pku.cn/) - A Shiny web-application accompanying https://doi.org/10.1038/s41597-019-0131-5. Source code available at https://github.com/Japrin/STARTRAC
- [Loupe Cell Browser](https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser), from 10x Genomics
- [Lung Aging Atlas](http://146.107.176.18:3838/MLAA_backup), reachable from https://theislab.github.io/LungAgingAtlas/ and linked to this publication: https://www.nature.com/articles/s41467-019-08831-9
- [Malaria Cell Atlas](https://www.sanger.ac.uk/science/tools/mca/mca/), to explore single parasite transcriptomes. Publication here: https://science.sciencemag.org/content/365/6455/eaaw2619
- [Mammary Gland Development](https://marionilab.cruk.cam.ac.uk/mammaryGland/), companion to https://www.nature.com/articles/s41467-017-02001-5
- [MERmaid](https://jef.works/MERmaid/), static HTML and Javascript browser for multiplexed FISH (MERFISH) data.
- [MicrobiomeDB](https://microbiomedb.org/mbio/app/), for exploring a large number of microbiome datasets. Original publication: https://academic.oup.com/nar/article/46/D1/D684/4584629
- [Microglia Orthologoues Viewer](https://amitlab.shinyapps.io/Orthologous_viewer/), for having easy data exploration and access to the resource presented in https://www.cell.com/cell/fulltext/S0092-8674(19)31231-0
- [Microglia Single Cell Atlas](http://www.microgliasinglecell.com), from the Stevens Lab. Related to manuscript https://www.cell.com/immunity/fulltext/S1074-7613(18)30485-0
- [Morpheus](https://software.broadinstitute.org/morpheus/), from the Broad Institute
- [Mousebrain.org](http://mousebrain.org), an atlas of cell types from the Linnarsson Lab - Manuscript: https://www.sciencedirect.com/science/article/pii/S009286741830789X
- [Mouse Cell Atlas](http://bis.zju.edu.cn/MCA/) Microwell-seq, companion to https://doi.org/10.1016/j.cell.2018.02.001
- [Mouse Gastrulation and Early Organogenesis](https://marionilab.cruk.cam.ac.uk/MouseGastrulation2018/), accompanying https://www.nature.com/articles/s41586-019-0933-9
- [MouseLight Neuron Browser](http://ml-neuronbrowser.janelia.org), an interactive web platform to explore, search, filter and visualize single neuron reconstructions generated by the Janelia MouseLight project. Related to https://www.biorxiv.org/content/10.1101/537233v1
- [Mouse Organogenesis](https://marionilab.cruk.cam.ac.uk/organogenesis/), companion to https://www.nature.com/articles/s41556-017-0013-z
- [Mouse Organogenesis Cell Atlas (MOCA)](https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/landing), companion to the manuscript https://www.nature.com/articles/s41586-019-0969-x
- [Mouse sci-ATAC-seq Atlas](http://atlas.gs.washington.edu/mouse-atac/), for example http://atlas.gs.washington.edu/mouse-atac/ukbb_visualization/, among others. Manuscript: https://doi.org/10.1016/j.cell.2018.06.052
- [Mouse Visual Cortex Explorer](http://greenberg.hms.harvard.edu/project/gene-database/), related to the manuscript https://doi.org/10.1038/s41593-017-0029-5 (Hrvatin et al., 2018)
- [Multiple sclerosis gene atlas](http://msatlas.compbio.sdu.dk), companion to https://www.biorxiv.org/content/10.1101/584920v1
- [NCI Patient Derived Model Repository](https://pdmdb.cancer.gov/pls/apex/f?p=101:1:0::NO) expression and whole exome data on >350 patients xenografts and patient tissue samples.
- [Neuroexpresso](http://neuroexpresso.org), for the analysis of brain cell type expression profiles. Source available here, https://github.com/PavlidisLab/neuroexpresso, publication: https://www.eneuro.org/content/4/6/ENEURO.0212-17.2017
- [nichExplorer](https://compbio.nyumc.org/niche/), an interactive browser for bone marrow microenvironment, accompanying https://doi.org/10.1038/s41586-019-1104-8, with source code available at https://github.com/igordot/nichexplorer. Also available from http://aifantislab.com/niche
- [OpenTargets Browser](https://www.opentargets.org/projects/effectorness), with datasets and apps on the effectorness gradient shaping the response of CD4+ T cells to cytokines. Publication available at https://www.nature.com/articles/s41467-020-15543-y, apps are https://cytokines.cellgeni.sanger.ac.uk/ and on cellxgene: https://cytokines.cellgeni.sanger.ac.uk/resting + https://cytokines.cellgeni.sanger.ac.uk/simulated
- [PanglaoDB](https://panglaodb.se), a database for the exploration of single cell RNA sequencing experiments from mouse and human. Publication at https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041, GitHub repo: https://github.com/oscar-franzen/PanglaoDB (not the source code of the resource itself)
- [PAWER](https://biit.cs.ut.ee/pawer/), a software tool for analysing protein microarray data - described in publication: https://www.biorxiv.org/content/10.1101/692905v1
- [PEAC RNA-seq](https://peac.hpc.qmul.ac.uk), for exploring RNA-seq of early rheumatoid arthritis individuals and their clinical data. Linked to the publication https://doi.org/10.1016/j.celrep.2019.07.091
- [Planaria Single Cell Atlas](https://shiny.mdc-berlin.de/psca/), related to the publication available at http://science.sciencemag.org/content/360/6391/eaaq1723
- [Planarian Digiworm](https://digiworm.wi.mit.edu), companion to http://science.sciencemag.org/content/360/6391/eaaq1736.full
- [pRolocGUI](http://bioconductor.org/packages/release/bioc/html/pRolocGUI.html), for interactive visualisation of spatial proteomics data
- [ProteinPaint](https://proteinpaint.stjude.org/examples/), some general info are here: https://www.stjude.org/research/shared-resources/technology-licensing/technologies/proteinpaint-web-application-for-visualizing-genomic-data-sj-15-0021.html, original publication here: https://www.nature.com/articles/ng.3466. Used e.g. to deploy Mouse Retina Development data, https://proteinpaint.stjude.org/F/mm10/2019.mouse.retina.scRNA.html
- [ProTrack](http://ccrcc.cptac-data-view.org), an interactive multiomics data browser for proteogenomic projects - described in the preprint https://doi.org/10.1101/2020.02.05.935650. Source code is (should be) available at https://github.com/WangLab/ProTrack-ccrcc (currently not found?)
- [puma_gtex](https://kuijjer.shinyapps.io/puma_gtex/), a browser for exploring tissue-specific functions of miRNAs (as described in the PUMA preprint, https://www.biorxiv.org/content/10.1101/2019.12.18.874065v1)
- [Reptilian pallium](https://public.brain.mpg.de/shiny/apps/SingleCellTranscriptomics/) More exotic species: https://brain.mpg.de/research/laurent-department/software-techniques.html
- [RNA-Seq Viewer](https://github.com/NCBI-Hackathons/rnaseqview) described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5605993/, and based on the [ideogram.js](https://github.com/eweitz/ideogram) library
- [RNA_vs_protein](https://tsomics.shinyapps.io/RNA_vs_protein/), to explore the quantitative proteome map of the human body, related to the publication https://www.biorxiv.org/content/10.1101/797373v2
- [scClustViz](https://baderlab.github.io/scClustViz/) - focused on exploring clustering results, paper here: https://f1000research.com/articles/7-1522/v2
- [scfind](https://scfind.sanger.ac.uk) for fast searches of large collections of single cell data - source code available here: https://github.com/hemberg-lab/scfind, described in the publication https://www.biorxiv.org/content/10.1101/788596v1
- [SCPortalen](http://single-cell.clst.riken.jp), a single-cell centric database. Published at https://doi.org/10.1093/nar/gkx949
- [scRNASeqDB](https://bioinfo.uth.edu/scrnaseqdb), a database for RNA-Seq gene expression profiles in human single cells. Publication https://www.mdpi.com/2073-4425/8/12/368, available at https://bioinfo.uth.edu/scrnaseqdb/
- [scSVA](https://github.com/klarman-cell-observatory/scSVA) for Single-cell Scalable Visualization and Analytics, presented in this publication https://www.biorxiv.org/content/10.1101/512582v1 - some documentation presented as videos, e.g. https://www.youtube.com/watch?v=W_U9swh468U
- [SCelVis](https://github.com/bihealth/scelvis), exploratory tool for preprocessed scRNA-seq data. Preprint: https://www.biorxiv.org/content/10.1101/713008v1 
- [SCope](http://scope.aertslab.org/), general visualization tool for large-scale scRNA-seq datasets. Has some datasets preprocessed to showcase, allows own data (in the .loom format) to be uploaded. Source code: https://github.com/aertslab/SCope. Original publication where this is found: https://www.sciencedirect.com/science/article/pii/S0092867418307207
- [SCV (Single Cell Viewer)](https://github.com/neuhausi/single-cell-viewer), preprint at https://www.biorxiv.org/content/10.1101/664789v2
- [Seattle Organismal Molecular Atlases](http://atlas.gs.washington.edu/hub/), the SOMA Data Portal
- [SeuratWizard](https://github.com/nasqar/SeuratWizard) - web-based interactive application to perform a guided single-cell RNA-seq data analysis and clustering based on Seurat. Based on the Nasqar platform (http://nasqar.abudhabi.nyu.edu)
- [SeuratV3Wizard](https://github.com/nasqar/seuratv3wizard) - successor of `SeuratWizard` (above) for Seurat v3
- [shinyMethyl](https://www.bioconductor.org/packages/release/bioc/html/shinyMethyl.html), Interactive visualization for Illumina methylation arrays
- [Single Cell Explorer](http://www.singlecellexplorer.org), published in https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6053-y. Demos available at http://54.159.6.229:8000/, source code at https://github.com/d-feng/scExplorer
- [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home), where the June 2019 release contains 99 single cell RNA-Seq studies (total: 882958 cells), across 11 different species. Publication: https://academic.oup.com/nar/article/47/D1/D711/5144130
- [Single Cell Portal](https://portals.broadinstitute.org/single_cell) at the Broad Institute; example here: https://portals.broadinstitute.org/single_cell/study/SCP43/atlas-of-human-blood-dendritic-cells-and-monocytes#study-visualize - Source code available here: https://github.com/broadinstitute/single_cell_portal_core
- [Slide-seq Portal](https://portals.broadinstitute.org/single_cell/study/SCP354/slide-seq-study), companion to https://science.sciencemag.org/content/363/6434/1463. Code available at https://github.com/broadchenf/Slideseq/tree/master/BeadSeq%20Code
- [Spaniel](https://bioconductor.org/packages/Spaniel/), designed to visualise results of Spatial Transcriptomics experiments; publication here https://www.biorxiv.org/content/10.1101/619197v1
- [SpatialDB](https://www.spatialomics.org/SpatialDB/), a database for spatially resolved transcriptomes. Described in the publication available at https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz934/5622703
- [spatialLIBD](http://research.libd.org/spatialLIBD/), a package and an app (http://spatial.libd.org/spatialLIBD/) to browse the LIBD human dorsolateral pre-frontal cortex (DLPFC) spatial transcriptomics data generated with the 10x Genomics Visium platform - described in the preprint https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1. Source available at https://github.com/LieberInstitute/spatialLIBD
- [ST Viewer](https://github.com/SpatialTranscriptomicsResearch/st_viewer), a GUI tool for easy and smooth visualisation and analysis of Spatial Transcriptomics datasets. Publication at https://academic.oup.com/bioinformatics/article/35/6/1058/5078471
- [TCGABrowser](http://tcgabrowser.ethz.ch:3838/PROD/), with source code available at https://github.com/pcheng84/TCGAbrowser/
- [TILAtlas](http://TILAtlas.org), the transcriptomic landscape of tumor-infiltrating CD8 lymphocytes in B16 melanoma tumors with single-cell RNA-Seq. Based on iSEE. Publication available at https://www.tandfonline.com/doi/full/10.1080/2162402X.2020.1737369 - alternative address: https://tilatlas.shinyapps.io/B16_CD8TIL_10X/
- [Transcriptome Dynamics for CD4+ T-cells in Experimental Malaria](https://malariaimmunology.shinyapps.io/cd4_dynamics/), developed as a companion to https://www.biorxiv.org/content/10.1101/675967v1
- [UCSC Cell Browser](https://cells.ucsc.edu/), a static HTML and Javascript browser for scRNA-seq data. Documentation at https://cellbrowser.readthedocs.io, described also in https://www.ncbi.nlm.nih.gov/pubmed/30407534
- [VDJView](https://bitbucket.org/kirbyvisp/vdjview), for analyzing immune single cell multi-omics data. Preprint available at https://www.biorxiv.org/content/10.1101/613083v2.full
- [Vidjil](http://www.vidjil.org/), for the visualization of High-Throughput Analysis
of V(D)J Immune Repertoire. Demo app here http://app.vidjil.org/, source code of the releases in http://www.vidjil.org/releases/
- [VisCello](https://cello.shinyapps.io/celegans/) for C. elegans embryogenesis, accompanying [https://www.biorxiv.org/content/10.1101/565549v2](https://www.biorxiv.org/content/10.1101/565549v2), with source code available at [https://github.com/qinzhu/VisCello](https://github.com/qinzhu/VisCello)
- [VirtualCytometry](https://www.grnpedia.org/cytometry/) from https://netbiolab.org/, focused on immune cell differentiation
- [Vitessce](http://vitessce.io), a visual integration tool for exploration of (spatial) single-cell experiment data. Some demos are available, e.g. http://vitessce.io/?dataset=linnarsson-2018 - source code is available at https://github.com/hubmapconsortium/vitessce
- [WIlsON](https://cran.r-project.org/web/packages/wilson/), described in https://academic.oup.com/bioinformatics/article/35/6/1055/5078467; demo here http://loosolab.mpi-bn.mpg.de
- [Xena Browser](https://xenabrowser.net), the UCSC Xena Functional Genomics Explorer
- [Zebrafish single cell development atlas](https://kleintools.hms.harvard.edu/paper_websites/wagner_zebrafish_timecourse2018/mainpage.html), related to the publication in http://science.sciencemag.org/content/360/6392/981

## Other resources:

- A collection of tools related to the Human Cell Atlas project: https://prod.data.humancellatlas.org/analyze
- An overview/comment on some of the above: http://science.sciencemag.org/content/358/6360/172
- A number of resources can be found at http://www.teichlab.org/data, although the newer sets adopted the use of iSEE

## Contributors

- Federico Marini [@federicomarini](https://github.com/federicomarini)
- Charlotte Soneson [@csoneson](https://github.com/csoneson)
- Kevin Rue-Albrecht [@kevinrue](https://github.com/kevinrue)
- Kevin Blighe [@kevinblighe](https://github.com/kevinblighe)
- Ben Busby [@DCGenomics](https://github.com/DCGenomics)
- David McGaughey [@davemcg](https://github.com/davemcg)
- Kamil Slowikowski [@slowkow](https://github.com/slowkow)
- Rajesh Patidar [@patidarr](https://github.com/patidarr)
- Ming Tang [@crazyhottommy](https://github.com/crazyhottommy)
- Michael von Papen [@mvonpapen](https://github.com/mvonpapen)
- Igor Dolgalev [@igordot](https://github.com/igordot)



