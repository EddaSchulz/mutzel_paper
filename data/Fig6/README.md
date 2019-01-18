# Structure of the data files in used in Fig6

1. *pyro\_tsix\_matlab.txt*:
Pyrosequencing data on XX. The file is structured as follows:
% 129 upstream | % 129 downstream
	- R1: Rep1 0h
	- R2: Rep1 0.5h
	- R3: Rep1 1h
	- R4: Rep1 2h
	- R5: Rep1 4h
	- R6: Rep1 8h
	- R7: Rep2 0h
	- R8: Rep2 0.5h
	- R9: Rep2 1h
	- R10: Rep2 2h
	- R11: Rep2 4h
	- R12: Rep2 8h
	- R13: Rep3 0h
	- R14: Rep3 0.5h
	- R15: Rep3 1h
	- R16: Rep3 2h
	- R17: Rep3 4h
	- R18: Rep3 8h

2. *130703\_matlab.txt*
qPCR data on doxycycline induced cells. The file is structured as follows:

	|       |      |Arpo | Gapdh | rrm2 | spliced Tsix | Tsix 3 | Xist | Tsix 5 | Rnf12 | Jpx | Tsx| 
	|------ |:-----|:-----|:-----:|:----:|:-----------: |:------:|:----:|:------:|:-----:|:----:|:----:|
	| **R1-12** | **Rep1 XO** | | 
	|  | R1-2: 0h |
	|  | R3-4: 0.5h |
	|  | R5-6: 1h |
	|  | R7-8: 2h |
	|  | R9-10: 4h |
	|  | R11-12: 8h |
	| **R13-24** | **Rep1 XX** |
	|  | R13-14: 0h |
	|  | R15-16: 0.5h |
	|  | R17-18: 1h |
	|  | R19-20: 2h |
	|  | R21-22: 4h |
	|  | R23-24: 8h |
	| **R25-36** | **Rep2 XO** |
	|  | R25-26: 0h |
	|  | R27-28: 0.5h |
	|  | R29-30: 1h |
	|  | R31-32: 2h |
	|  | R33-34: 4h |
	|  | R35-36: 8h |
	| **R37-48** | **Rep2 XX** |
	|  | R37-38: 0h |
	|  | R39-40: 0.5h |
	|  | R41-42: 1h |
	|  | R43-44: 2h |
	|  | R45-46: 4h |
	|  | R47-48: 8h |
	| **R49-60** | **Rep3 XO** |
	|  | R49-50: 0h |
	|  | R51-52: 0.5h |
	|  | R53-54: 1h |
	|  | R55-56: 2h |
	|  | R57-58: 4h |
	|  | R59-60: 8h |
	| **R61-72** | **Rep3 XX** |
	|  | R61-62: 0h |
	|  | R63-64: 0.5h |
	|  | R65-66: 1h |
	|  | R67-68: 2h |
	|  | R69-70: 4h |
	|  | R71-72: 8h |
	
	2 subsequent rows are duplicates for the same sample.	

3. *quant\_TXY\_Rep1\_all.txt*, *quant\_TXY\_Rep1\_bgd.txt*
Results of the RNA-FISH signal quantification with Xist, Tsix3\' and DxPas34 probes.
*_all* contains the quantification of true transcription sites while *_bgd* contains the quantification of random nuclear background ROIs.
Analogue for *Rep2* and *Rep3*. The files are structured as follows:

	|C1: Channel | C2: ROI quantification| 
	|------------ | -------------|
	Three consecutive rows quantify the same ROI in the three different channels.\
	Channel 1: Xist\
	Channel 2: DxPas34\
	Channel 3: Tsix3´
	
4. The folder *Microscopy_images* contains the Microscopy images of Fig6.
