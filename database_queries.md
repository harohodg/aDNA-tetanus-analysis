# Query to search for all runs with tetanus identified in the STAT analysis
The EXCEPT clause removes columns with data that can't be exported to csv

```
SELECT m.* 
    EXCEPT 
    (attributes, 
    biosamplemodel_sam, 
    geo_loc_name_sam, 
    ena_first_public_run, 
    ena_last_update_run, 
    sample_name_sam, 
    jattr,
    datastore_filetype,
    datastore_provider,
    datastore_region), 
    tax.total_count, 
    tax.self_count 
    FROM nih-sra-datastore.sra.metadata as m, 
    nih-sra-datastore.sra_tax_analysis_tool.tax_analysis as tax 
    WHERE m.acc=tax.acc and tax_id=1513
    ORDER BY tax.total_count
```

Could be simplified to 
```
SELECT m.acc, tax.self_count, tax.total_count 
    FROM nih-sra-datastore.sra.metadata as m, 
    nih-sra-datastore.sra_tax_analysis_tool.tax_analysis as tax 
    WHERE m.acc=tax.acc and tax_id=1513
    ORDER BY tax.total_count DESCENDING
```
and then use the [SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/) to obtain the corresponding meta data


# Query to get summary of SRA database
```
SELECT 
COUNT(DISTINCT acc) as runs, 
COUNT(DISTINCT experiment) as experiments,
COUNT(DISTINCT biosample) as biosamples,
COUNT(DISTINCT sra_study) as studies,
COUNT(DISTINCT bioproject) as bioprojects,
SUM( mbytes ) as totalMegaBytes,
SUM( mbases ) as totalMegaBases,
CURRENT_DATE() as the_date
FROM `nih-sra-datastore.sra.metadata`
```

# Query to get full STAT analysis for each run
```
SELECT  * FROM
    nih-sra-datastore.sra_tax_analysis_tool.tax_analysis 
    AS tax 
    WHERE tax.acc 
    IN
    ("DRR046399",
    "DRR046400",
    "DRR046402",
    "DRR046405",
    "DRR046409",
    "DRR046398",
    "DRR046401",
    "DRR046403",
    "DRR046404",
    "DRR046408",
    "SRR1238558",
    "SRR1238559",
    "SRR9898505",
    "SRR5047069",
    "SRR1314212",
    "SRR1314214",
    "SRR1298752",
    "ERR4374869",
    "ERR4374914",
    "ERR4375049",
    "ERR4375137",
    "ERR4375141",
    "ERR2862150",
    "ERR2862151",
    "ERR3828670",
    "ERR3841646",
    "ERR966294",
    "ERR966295",
    "ERR966296",
    "ERR966297",
    "ERR966298",
    "ERR966302",
    "ERR966303",
    "ERR966304",
    "ERR966305",
    "ERR966306",
    "ERR966307",
    "ERR966312",
    "ERR966313",
    "ERR966426",
    "SRR5169887",
    "SRR5169888",
    "ERR3426230",
    "ERR1193532",
    "ERR651004",
    "ERR1368878",
    "ERR2111944",
    "ERR2111945",
    "ERR2111946",
    "ERR2111947",
    "ERR2111948",
    "ERR2111949",
    "ERR2111950",
    "ERR2111951",
    "ERR2112074",
    "ERR2112075",
    "ERR2112076",
    "ERR2112077",
    "ERR2112078",
    "ERR2112079",
    "ERR2112080",
    "ERR2112081",
    "ERR2112091",
    "ERR2112098",
    "ERR2225784",
    "ERR2204614",
    "ERR2112574",
    "ERR2112575",
    "ERR2112576",
    "ERR2112577",
    "ERR1937833",
    "ERR1937835",
    "ERR2270785",
    "ERR4017942",
    "ERR4017944",
    "ERR4017946")
```
