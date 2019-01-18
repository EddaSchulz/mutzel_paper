# Structure of the data files in used in Fig2,3,4

1. *Data of E5.0.xlsx*:\
Xist/Tsix RNA FISH measurements in female mouse epiblast cells at E5.0 of embryogenesis. The file is structured as follows:

	|Embryo # | cells with 2xXist clouds | cells with 1xXist cloud  | cells with 0xXist clouds|
	|--|---|---|---|

	Each row is one embryo.

1. *biallelic_FISH_exp1.txt*, *biallelic_FISH_exp2.txt*, *biallelic_FISH_exp9.txt*:\ 
3 replicates of the artificial biallelic induction experiment (Xist RNA-FISH measurements). Each file is structured as follows:

| | 48h prior to Dox | 0h prior to Dox | 6h -Dox | 24h -Dox | 48h -Dox | 6h +Dox | 24h +Dox | 48h +Dox|
|------------ | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------|
XX cells 2xXist |  |  |  |  |  |  |  |  | 
XX cells 1xXist |  |  |  |  |  |  |  |  | 
XX cells 0xXist |  |  |  |  |  |  |  |  | 
XO cells 1xXist |  |  |  |  |  |  |  |  | 
XO cells 0xXist |  |  |  |  |  |  |  |  | 
cells w/o any signal|  |  |  |  |  |  |  |  | 

1. Published data\
	The data have the following structure:\
	
	|time point | cells w/o Xist signal | cells with 1 Xist signal | cells with 2 Xist signals|
	|--|---|---|---|
	
	1. *mouse_in_vivo_data.txt*: Sakata et al, Development 144: 2784–2797 (2017) and this study
	1. *rabbit_in_vivo_data_pp.txt*: Okamoto et al, Nature 472: 370–374 (2011)
	
1. *TX_fish_a.txt*, *TX_fish_b.txt*, *TX_fish_c.txt*\
3 Replicates of the doxycycline speedup experiment (Xist RNA-FISH measurements). Each file is structured as follows:

	| | 2xXist | 1xXist | 0xXist XX | XO| 
	|------------ | ------------ | ------------ | ------------ | ------------|
	0d of Diff + Dox |  |  |  | 
	1d of Diff + Dox |  |  |  | 
	2d of Diff + Dox |  |  |  | 
	3d of Diff + Dox |  |  |  | 
	4d of Diff + Dox |  |  |  | 
	0d of Diff -Dox |  |  |  | 
	1d of Diff -Dox |  |  |  | 
	2d of Diff -Dox |  |  |  | 
	3d of Diff -Dox |  |  |  | 
	4d of Diff -Dox |  |  |  | 

1. *xist_fold.txt*\
Illumina targeted expression assay (TRex). The file is structured as follows:

	time point | replicate | b6 or cast | Xist Dox/Cntrl 
	|--|---|---|---|

	Each Row is 1 time point in one replicate

1. *results_96h_IF_FISH.txt*:\
H3K27me3 immunofluorescence and Xist RNA-FISH 48h after doxycycline induction (96h of differentiation). The file is structured as follows:

	| | 96h -Dox Rep1 | 96h +Dox Rep1 | 96h -Dox Rep2 | 96h +Dox Rep2 | 96h -Dox Rep3 | 96h +Dox Rep3 |
	|------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------|
	|% XiXa |  |  |  |  |  | 
	|% XiXa\* |  |  |  |  |  | 
	|% XiXi |  |  |  |  |  | 
	|% divers |  |  |  |  |  |   

1. *edu_data_table.txt*, *edu_additional_data_table.txt*:\
The first file contains counts of randomly selected cells. To obtain more biallelically Xist expressing cells these were counted selectively and are saved in the second file.
For the analysis both counts were summed up.
Edu Labeling together with Xist RNA-FISH 48h after doxycycline induction (96h of differentiation). The files are structured as follows:

	| | 54h -Dox Rep1 | 54h +Dox Rep1 | 72 -Dox Rep1 | 72h +Dox Rep1 | 96h -Dox Rep1 | 96h +Dox Rep1 | 54h -Dox Rep2 | 54h +Dox Rep2 | 72h -Dox Rep2 | 72h +Dox Rep2 | 96h -Dox Rep2 | 96h +Dox Rep2 | 54h -Dox Rep3 | 54h+Dox Rep3 | 72h -Dox Rep3 | 72h +Dox Rep3 | 96h -Dox Rep3 | 96h +Dox Rep3|
	|------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------
	|EdU+ cells 1xXist |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 
	|EdU+ cells 2xXist |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 
	|EdU+ cells 0xXist |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 
	|EdU- cells 1xXist |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 
	|EdU- cells 2xXist |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 
	|EdU- cells 0xXist |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 

1. The folder *Microscopy_images* contains the Microscopy images of Fig3, Fig4 and FigS3.
