import extractdata
import refgenome
import pipeline

test_genome = refgenome.ReferenceGenome("/Users/jbreynier/Desktop/Research/Yang Lab/Fall Quarter/hg38/hg38.fa",
                                        "/Users/jbreynier/Desktop/Research/Yang Lab/Fall Quarter/hg38/hg38.chrom.sizes",
                                        "/Users/jbreynier/Desktop/Research/Yang Lab/Fall Quarter/hg38/k100.umap.bed")
test_pipeline = pipeline.MotifPipeline()
test_pipeline.set_output_dir("/Users/jbreynier/Desktop/Research/Yang Lab/Fall Quarter/Yang_Yang_analysis/test_pipeline_2")
test_pipeline.set_input_dir("/Users/jbreynier/Desktop/Research/Yang Lab/Fall Quarter/Yang_Yang_analysis/PBTA_bedpe")
test_pipeline.set_rand_sv_ratio("3")
test_pipeline.set_SV_types("del,dup,inv,tra")
test_pipeline.list_bedpe = ["BS_ZPF89720", "BS_0DKPGQWD", "BS_0EN5R1AP"]

extractdata.bedpe_to_bed(test_genome, test_pipeline)