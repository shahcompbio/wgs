### WGS qc pipelines docs

#### WGS cohort qc

1. The format of the input yaml for cohort_qc is extremely simple:
```
cohort_name:
  patient_name:
    sample_maf: {MAF}
    sample_label: {LABEL}
    ...
  patient_name:
    sample_maf: {MAF}
    sample_label: {LABEL}
  ...
```
The sample_label field is what is used as the patient label in the pipeline's output, and can be used to add sample-level groupings to outputs.

The input {MAF} can be either the germline or consensus maf outputted from `wgs variant_calling`

2. The pipeline uses the cBioPortal API and requires an api key to run the pipeline.
``` 
wgs cohort_qc --input_yaml {input_yaml} --API_key {key} --outdir {out} --tmppdir {tmp}
```
The API key can be obtained from [this website](https://docs.cbioportal.org/6.-web-api-and-clients/api-and-api-clients).

3. In order to run the pipeline you also need to use R scripts installable in this [conda recipe](https://github.com/shahcompbio/conda-recipes/tree/master/pseudo_bulk_qc_html_report)
#### WGS sample qc
1. The format of the input yaml for wgs sample qc:
```
Patient:
  breakpoints_consensus: {filtered_consensus_calls file from wgs variant_calling 0.0.1 or 0.1.0}
  germline_calls:  {germline vcf or maf from wgs variant calling 0.0.1 or 0.1.0}
  normal_bam: {patient normal bam}
  roh: {roh from wgs variant_calling 0.0.1 or 0.1.0}
  somatic_calls:  {germline vcf or maf from wgs variant calling 0.0.1 or 0.1.0}
  titan: {titan markers from wgs copynumber 0.1.0 and 0.0.1 }
  tumour_bam: {patient tumor bam}
  remixt: {remixt h5 from wgs wgs copynumber 0.1.0 and 0.0.1}
```

2. the command format is

```
wgs cohort_qc --input_yaml {input_yaml} --outdir {out} --tmppdir {tmp}```