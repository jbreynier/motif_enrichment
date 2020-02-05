import sys
import getopt
import refgenome
import pipeline
import exceptions
import pipeline
import extractdata
import runprogram
import graphs

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                "i:o:f:l:e:m:a:t:s:r:F:A:h",
                                ["input_dir=",
                                    "output_dir=",
                                    "genome_fasta=",
                                    "genome_len=",
                                    "genome_include=",
                                    "motif_path=",
                                    "sample_attr=",
                                    "sampleinfo_table=",
                                    "SV_types=",
                                    "rand_sv_ratio=",
                                    "FIMO_thresh=",
                                    "AME_scoring=",
                                    "help"
                                    ]
                                )
    except getopt.GetoptError as e:
        print(e)
        sys.exit(2)
    
    if len(args) > 0:
        message = "Error: non-paired arguments are not allowed."
        raise exceptions.WrongArgumentError(message)
    
    motif_pipeline = pipeline.MotifPipeline()
    sample_attr_path = None
    genome_fasta = None
    genome_len = None
    genome_include = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            description()
            sys.exit()
        elif opt in ("-i", "--input_dir"):
            motif_pipeline.set_input_dir(arg)
        elif opt in ("-o", "--output_dir"):
            motif_pipeline.set_output_dir(arg)
        elif opt in ("-f", "--genome_fasta"):
            genome_fasta = arg
        elif opt in ("-l", "--genome_len"):
            genome_len = arg
        elif opt in ("-e", "--genome_include"):
            genome_include = arg
        elif opt in ("-m", "--motif_path"):
            motif_pipeline.set_motif_path(arg)
        elif opt in ("-a", "--sample_attr"):
            motif_pipeline.set_sample_attr(arg)
        elif opt in ("-t", "--sampleinfo_table"):
            sample_attr_path = arg
        elif opt in ("-s", "--SV_types"):
            motif_pipeline.set_SV_types(arg)
        elif opt in ("-r", "--rand_sv_ratio"):
            motif_pipeline.set_rand_sv_ratio(arg)
        elif opt in ("-F", "--FIMO_thresh"):
            motif_pipeline.set_FIMO_thresh(arg)
        elif opt in ("-A", "--AME_scoring"):
            motif_pipeline.set_AME_scoring(arg)
        else:
            message = "Error: {opt} is not a valid option".format(opt=opt)
            raise exceptions.WrongArgumentError(message)
    
    if sample_attr_path is None:
        message = "Error: you must indicate --sampleinfo_table."
        raise exceptions.MissingArgumentError(message)
    if genome_fasta is None:
        message = "Error: you must indicate --genome_fasta."
        raise exceptions.MissingArgumentError(message)
    if genome_len is None:
        message = "Error: you must indicate --genome_len."
        raise exceptions.MissingArgumentError(message)
    if genome_include is None:
        message = "Error: you must indicate --genome_include."
        raise exceptions.MissingArgumentError(message)
    for pipeline_attr in ["input_dir",
                            "output_dir",
                            "motif_path",
                            "sample_attr",
                            "SV_types",
                            "rand_sv_ratio",
                            "FIMO_thresh",
                            "AME_scoring"]:
        if not hasattr(motif_pipeline, pipeline_attr):
            message = ("Error: you must indicate --{attr}.").format(attr=pipeline_attr)
            raise exceptions.MissingArgumentError(message)
    motif_pipeline.write_description()
    motif_pipeline.set_list_bedpe(sample_attr_path)
    reference_genome = refgenome.ReferenceGenome(genome_fasta, genome_len, genome_include)
    extractdata.bedpe_to_bed(reference_genome, motif_pipeline)
    motif_pipeline.set_num_SV_breakpoints()
    runprogram.bedtools(motif_pipeline, reference_genome)
    runprogram.FIMO(motif_pipeline)
    runprogram.AME(motif_pipeline)
    extractdata.extract_list_sequences_AME(motif_pipeline)
    extractdata.extract_output_FIMO(motif_pipeline)
    extractdata.extract_output_AME(motif_pipeline)
    graphs.generate_histogram(motif_pipeline)

def description():
    manual = ("\nOPTIONS:\n\t -h or --help : display the manual\n"
        "\t -d or --g_dir : specify the genome directory path\n"
        "\t -f or --g_fa: specify the hg19 genome fasta file path\n"
        "\t -g or --g_gtf : specify the hg19 genome GTF file path\n"
        "\t -r or --r1r2_dir : specify the fastq directory for the sample\n"
        "\t -a or --align_dir : specify the directory of "
        "the aligned bam file\n"
        "\t -o or --out : specify the final output directory for TRUST\n"
        "\t -p or --param : specify additional parameters for TRUST\n"
        "\t -s or --log_suffix : specify suffix for log file name\n"  
        "\t -l or --loop : input for -r/--r1r2_dir and -a/--align_dir \n"
        "\t are directories containing different samples to loop over\n"
        )
    print(manual)

if __name__ == "__main__":
    main()