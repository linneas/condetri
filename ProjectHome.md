During the last few years, DNA and RNA sequencing have become to play a more important role in biological and medical applications, especially due to the increasing amount of sequencing data yield from the sequencing machines and the enormous drop down of sequencing costs. Particularly, Illumina sequencing has an increasing impact on gathering data from model and non-model organisms. However, accurate and easy to use tools for quality filtering have not yet been established.

We present ConDeTri, a method for content dependent read trimming for Illumina data using quality scores of each base individually. It is independent from sequencing coverage and user interaction. The main focus of the implementation is on usability and to incorporate read trimming in next-generation sequencing data processing and analysis pipelines. It can process single-end and paired-end sequencing data of arbitrary length.

ConDeTri is able to trim and remove reads with low quality scores to get better results from de novo assemblies. Especially, low coverage or large genome sequencing projects gain from trimming reads.

See the wiki page for more information and a description of the program.

Citation:<br>
Smeds L, Künstner A (2011) ConDeTri - A Content Dependent Read Trimmer for Illumina Data. PLoS ONE 6(10): e26314. doi:10.1371/journal.pone.0026314