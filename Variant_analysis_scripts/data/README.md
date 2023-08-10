To reproduce the variant analysis, you will need to download the read alignment tarball and extract the contents here.

You can get the tarball from:
[https://doi.org/10.6084/m9.figshare.23925768](https://doi.org/10.6084/m9.figshare.23925768)

Inspect its hash and make sure the downloaded archive matches:
`cat read_alignments.tar.gz.md5`
`930575976015ec307866bba2a8c1943a  read_alignments.tar.gz`

`md5sum read_alignments.tar.gz`
`930575976015ec307866bba2a8c1943a  read_alignments.tar.gz`

Decompress the archive with:
`tar -zxvf read_alignments.tar.gz`

The extracted contents should contain:
* `./read_alignments`, which contains:
  * `./read_alignments/chromosome` - chromosomal read alignments, in BAM format
  * `./read_alignments/plasmid` - plasmid read alignments
  * NOTE that both of these folders should contain sub-folders named `merged` and `singles`. I.e., there must be `./read_alignments/chromosome/merged`, `./read_alignments/chromosome/singles`, `./read_alignments/plasmid/merged`, and `./read_alignments/plasmid/singles`.

Also, decompress the reference genome archive:
`tar -zxvf ref_genome.tar.gz`

Then, you should be able to reproduce the analysis using: `../scripts/variants.sh --variants --alignments --threads [number of threads]`
