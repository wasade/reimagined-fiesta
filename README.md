# reimagined-fiesta
Expand deblur to closed reference

This code was originally put together to handle the situation of needing to run output from [deblur](https://github.com/biocore/deblur) against [PICRUSt](http://picrust.github.io/picrust/). 

The code takes as input, an OTU map relating the deblurred sequences to Greengenes the 13\_8 99% otus (i.e., what the PICRUSt reference is based off) and the deblur OTU table. The resulting table is in the namespace of Greengenes, and the values are the sum of the counts of a deblurred sequence in each sample for each deblurred sequence within a GGID in that sample. Or, a GG ID may contain multiple deblurred sequences, so we aggregate the original deblurred counts. 

```
python expand.py deblurred_input.biom otu_map.txt expanded_output.biom
```
