STAT 540 Project: Expression Profiling in Inflammatory Bowel Disease
====================

### Group Members and Specialties

<table><thead>
<tr>
<th>Name</th>
<th>Program</th>
<th>Supervisor</th>
<th>Expertise</th>
</tr>
</thead><tbody>
<tr>
<td>Yue Sun</td>
<td>Bioinformatics(PhD)</td>
<td>Jane Z. Wang</td>
<td>Bioinformatics</td>
</tr>
<tr>
<td>Yiming Zhang</td>
<td>Electrical and Computer Engineering(MSc)</td>
<td>Jane Z. Wang</td>
<td>Computer Science</td>
</tr>
<tr>
<td>Omar AlOmeir</td>
<td>Computer Science (MSc)</td>
<td>Rachel Pottinger</td>
<td>Data Management</td>
</tr>
<tr>
<td>Abrar Wafa</td>
<td>Electrical and Computer Engineering (MSc)</td>
<td>Rabab Ward</td>
<td>Signal Processing</td>
</tr>
</tbody></table>

### Definitions and Abbreviations

IBD: Inflammatory bowel disease, CD: Crohn's disease, UC: ulcerative colitis.

### Background

<!-- couple sentences of biological/scientific context
motivate interest in a broad line of inquiry -->

Inflammatory bowel disease (IBD) is a group of inflammatory conditions of the colon and small intestineis. And, it is a complex disease which arises as a result of the interaction of environmental and genetic factors. There are two main forms of IBD: Crohn's disease (CD) and ulcerative colitis (UC). There is an apparent overlap between the two forms in clinical criteria, understanding of pathophysiology, and therapy.  The study from which this data set originated [2] aimed to broaden the understanding of gene regulation events in CD and UC at a genome-wide level, and identify novel unknown genes involved in perpetuating inflammatory disease progression. 

### Data Description

<!-- 
0. Source of data?
1. what is being measured?
2. with what platform?
3. how many samples? 
4. are replicates biological or technical? How many replicates are there?-->

0. Public data set (GEO accession number：GSE1710)[1].
1. Differentially regulated genes.
2. High-density cDNA microarrays?
3. Platform :	GPL284, Human UniGene Set RZPD 1
4. 31 samples, split into 3 groups: 11 control samples, 10 Crohn's disease patients, 10 ulcerative colitis patients. Each sample is from a different individual (patient or control).
5. Biological? 34561 replicates.

### Work Plan

<!-- outline of analyses you can probably do with this data to answer those questions
GET SPECIFIC, at least as specific as you can
"identify differentially expressed genes between the wild type and the knockout" is better than "conduct appropriate statistical analyses"
even better is to enhance with more specifics, e.g. you anticipate using a linear model as implemented in limma -->

<table><thead>
<tr>
<th>Work item</th>
<th>Responsible team member(s)</th>
</tr>
</thead><tbody>
<tr>
<td>Data Intake. More specifically: loading, cleaning, and sorting the data.</td>
<td>Abrar Wafa</td>
</tr>
<tr>
<td>Quality Control. Exploring the data, looking for batch effects, outliers, et cetera.</td>
<td>Omar AlOmeir and Yiming Zhang</td>
</tr>
<tr>
<td>Differential Expression Analysis to detect differentially regulated genes. Namely: normal controls vs at least one of the IBD subtypes, normal vs CD, normal vs UC</td>
<td>Omar AlOmeir and Abrar Wafa</td>
</tr>
<tr>
<td>PCA and Clustering. Find a compact representation of the data, correct potential batch effects[3].</td>
<td>Yiming Zhang</td>
</tr>
<tr>
<td>Classification and Cross Validation. Classify genes in different groups, define features in expressions among samples, test and improve the classifier using cross validation.</td>
<td>Yue Sun and Yiming Zhang</td>
</tr>
</tbody></table>
<td>Gene Set Enrichment Analysis. Assign to functional groups based on classification by Gene Ontology </td>
<td>Yue Sun</td>
</tr>
<tr>

### References

[1] [Expression profiling in inflammatory bowel disease data set](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1710)

[2] [Costello, Christine M., et al. "Dissection of the inflammatory bowel disease transcriptome using genome-wide cDNA microarrays." PLoS medicine 2.8 (2005): e199.](http://www.plosmedicine.org/article/info%3Adoi%2F10.1371%2Fjournal.pmed.0020199#pmed-0020199-g004)

[3] [Chen, Chao, et al. "Removing batch effects in analysis of expression microarray data: an evaluation of six batch adjustment methods." PloS one 6.2 (2011): e17238.](http://www.plosone.org/article/fetchObject.action?uri=info%3Adoi%2F10.1371%2Fjournal.pone.0017238&representation=PDF)
