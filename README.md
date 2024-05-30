# Attentional Gait Index

The code used to generate the attentional gait index and reproduce results shown in Liu, S., Rosso, A. L., Baillargeon, E. M., Weinstein, A. M., Rosano, C., & Torres-Oviedo, G. (2024). Novel attentional gait index reveals a cognitive ability-related decline in gait automaticity during dual-task walking. Frontiers in aging neuroscience, 15, 1283376.https://www.frontiersin.org/articles/10.3389/fnagi.2023.1283376/full

## How do I create the attentional gait index for my dual-task study with fNIRS? 
Call the function computeAutomaticityIndex.m
The only required arguments are the perfDataPath, pfcDataPath, and studyID. The two paths of where the PFC data and performance data are stored. It requires the data to be follow a specific column format and header.
Note that the index requires at least two levels of difficulty.

## To reproduce the results in the figure
Run AutoIndexMainAnalysis.m
This file currently contains everything that's shown in the paper, and more.

> [!NOTE]
>Minor improvements on docs and organization is incoming. The codebase includes everything that's needed to reproduce the published results, and a lot more. Improvements could be made in the future to provide better docs (particularly on what file formatting are needed to recreate the index for other task designs), and separate exploratory analysis from the primary analysis included in the article. 
