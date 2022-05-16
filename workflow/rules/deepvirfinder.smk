##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule deepvirfinder:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "{sample}.filtered.sorted.fasta",
        ),
    output:
        txt=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "deepvirfinder",
            "{sample}",
            "{sample}_gt{cutoff}bp_dvfpred.txt",
        ),
    params:
        output_dir=directory(
            os.path.join(
                OUTPUT_FOLDER,
                "processing_files",
                "deepvirfinder",
                "{sample}",
            )
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "deepvirfinder",
            "{sample}.{cutoff}.deepvirfinder.log",
        ),
    resources:
        cpus=5,
    conda:
        "../envs/deepvirfinder.yaml"
    threads: 5
    shell:
        """
        if [[ "$OSTYPE" == "darwin"* ]]; then
            wget https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX10.9.sdk.tar.xz &> "{log}"
            tar -xf MacOSX10.9.sdk.tar.xz &> "{log}"
            export CONDA_BUILD_SYSROOT=MacOSX10.9.sdk
        fi

        dvf.py -i {input.contig} -o {params.output_dir} \
        -l {wildcards.cutoff} -c {threads} &> "{log}"

        mv "{params.output_dir}"/*_gt*txt "{output.txt}"

        if [[ "$OSTYPE" == "darwin"* ]]; then
            rm -rf MacOSX10.9*
        fi
        """


##########################################################################
##########################################################################
