import sys
import os
import glob
import getopt
import refgenome
import pipeline
import exceptions
import pipeline
import extractdata
import runprogram
import importlib
import subprocess
import graphs
import parsl

def main():
    '''Reads input from terminal and coordinates pipeline'''
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                "i:o:f:l:e:m:a:t:s:r:F:A:c:p:h",
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
                                    "config=",
                                    "prefix=",
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
    prefix = None
    config_name = "local"

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
        elif opt in ("-c", "--config"):
            config_name = arg
        elif opt in ("-p", "--prefix"):
            prefix = arg
        else:
            message = "Error: {opt} is not a valid option".format(opt=opt)
            raise exceptions.WrongArgumentError(message)
    
    if ((sample_attr_path is None and not motif_pipeline.sample_attr == "all") or 
        (sample_attr_path is not None and motif_pipeline.sample_attr == "all")):
        message = "Error: you must indicate both --sampleinfo_table and --sample_attr, or neither."
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
                            "motif_path"]:
        if not hasattr(motif_pipeline, pipeline_attr):
            message = ("Error: you must indicate --{attr}.").format(attr=pipeline_attr)
            raise exceptions.MissingArgumentError(message)
    motif_pipeline.set_subdir_name(prefix)
    motif_pipeline.write_description()
    motif_pipeline.set_list_bedpe(sample_attr_path)
    reference_genome = refgenome.ReferenceGenome(genome_fasta, genome_len, genome_include)
    base_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
    try:
        config = os.path.join(base_dir, 'configs', '{}.py'.format(config_name))
        spec = importlib.util.spec_from_file_location('', config)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        parsl.load(module.config)
    except:
        raise exceptions.IncorrectPathError("Cannot find the config file <{config_name}>.".format(config_name=config_name))
    
    if not os.path.isdir(motif_pipeline.output_dir+"bed_files"):
        os.mkdir(motif_pipeline.output_dir+"bed_files")
    for file_name in motif_pipeline.list_bedpe:
        sv_types_to_run = get_SV_types(motif_pipeline, file_name)
        if sv_types_to_run:
            extractdata.bedpe_to_bed(reference_genome, motif_pipeline, file_name, sv_types_to_run)
    parsl.wait_for_current_tasks()
    runprogram.merge(motif_pipeline)
    motif_pipeline.set_num_SV_breakpoints()
    runprogram.bedtools(motif_pipeline, reference_genome)
    runprogram.FIMO(motif_pipeline)
    runprogram.AME(motif_pipeline)
    extractdata.extract_list_sequences_AME(motif_pipeline)
    extractdata.extract_output_FIMO(motif_pipeline)
    extractdata.extract_output_AME(motif_pipeline)
    graphs.generate_histogram(motif_pipeline)

def get_SV_types(motif_pipeline, sample_name):
    '''Determines which SV type to run the analysis for for each sample'''
    sv_types_to_run = []
    file_prefix = motif_pipeline.output_dir + "bed_files/" + sample_name
    for sv_type in motif_pipeline.SV_types:
        if (not os.path.isfile(file_prefix+"_"+sv_type+"_sv.bed") or 
            not os.path.isfile(file_prefix+"_"+sv_type+"_rand"+str(motif_pipeline.rand_sv_ratio)+".bed")):
            sv_types_to_run.append(sv_type)
    return sv_types_to_run

def description():
    '''Prints description of pipeline usage'''
    manual = ("\nOPTIONS:\n\t -h or --help : display the manual and exit\n"
        "\t -i or --input_dir : specify the directory path containing all the input bedpe files\n"
        "\t -o or --output_dir : specify the directory for all output files\n"
        "\t -f or --genome_fasta: specify the genome fasta file path\n"
        "\t -l or --genome_len : specify the chromosome length file path\n"
        "\t -e or --genome_include : specify the file path of the genome regions to include\n"
        "\t -m or --motif_path : specify the path of the MEME motif file\n"
        "\t -a or --sample_attr : specify the sample attribute(s) to select for "
        "[format: < attribute_type:specific_attribute,attribute_type:specific_attribute >]\n"
        "\t -t or --sampleinfo_table : specify the sample attribute table file path\n"
        "\t -s or --SV_types : specify SV type(s) [format: < tra,inv,del,dup >]\n"
        "\t -r or --rand_sv_ratio : specify the ratio of random SVs to real SVs [format: < int >]\n"
        "\t -F or --FIMO_thresh : specify the p-value threshold for the FIMO algorithm [format: < float >]\n"
        "\t -A or --AME_scoring : specify the scoring method for AME [either < avg > or < max >]\n"
        "\t -c or --config : specify the config file to use for Parsl parallel processing\n"
        "\t -p or --prefix : specify a custom prefix for the results directory and files\n"
        )
    print(manual)

if __name__ == "__main__":
    main()