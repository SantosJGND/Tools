class ConstantsSettings:
    """
    classdocs
    """

    project_directory = "/tmp/televir/projects/"
    media_directory = "/tmp/televir/media/"
    static_directory = "/tmp/televir/static/"
    televir_subdirectory = "televir_projects"
    PAGINATE_NUMBER = 10
    docker_app_directory = "/televir/mngs_benchmark/"
    docker_install_directory = "/televir/mngs_benchmark/mngs_environments/"
    USER_TREE_INDEX = 0

    PIPELINE_NAME_read_quality_analysis = "Read quality analysis and improvement"
    PIPELINE_NAME_type_and_subtype_analysis = "Classification"
    PIPELINE_NAME_variant_detection = "Mutation detection and consensus generation"
    PIPELINE_NAME_coverage_analysis = "Coverage analysis"
    PIPELINE_NAME_alignment = "Alignment/Phylogeny"
    PIPELINE_NAME_intra_host_minor_variant_detection = (
        "Intra-host minor variant detection"
    )
    PIPELINE_NAME_viral_enrichment = "Viral enrichment"
    PIPELINE_NAME_host_depletion = "Host depletion"
    PIPELINE_NAME_contig_classification = "Contig classification"
    PIPELINE_NAME_read_classification = "Read classification"
    PIPELINE_NAME_assembly = "Assembly"
    PIPELINE_NAME_remapping = "Remapping"

    ## values to upload to database
    vect_pipeline_names = [
        PIPELINE_NAME_read_quality_analysis,
        PIPELINE_NAME_type_and_subtype_analysis,
        PIPELINE_NAME_variant_detection,
        PIPELINE_NAME_coverage_analysis,
        PIPELINE_NAME_alignment,
        PIPELINE_NAME_intra_host_minor_variant_detection,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_contig_classification,
        PIPELINE_NAME_read_classification,
        PIPELINE_NAME_assembly,
        PIPELINE_NAME_remapping,
    ]

    ###############################
    ### technology available
    TECHNOLOGY_illumina_old = "Illumina"
    TECHNOLOGY_illumina = "Illumina/IonTorrent"
    TECHNOLOGY_minion = "ONT"

    ###################################

    PIPELINE_STEPS_DB_DEPENDENT = [
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_read_classification,
        PIPELINE_NAME_contig_classification,
    ]

    ################################### Project_file_structure

    DIRS = {
        PIPELINE_NAME_read_quality_analysis: "reads/clean/",
        "reads_depleted_dir": "reads/hd_filtered/",
        "reads_enriched_dir": "reads/enriched/",
        PIPELINE_NAME_host_depletion: "host_depletion/",
        PIPELINE_NAME_viral_enrichment: "enrichment/",
        PIPELINE_NAME_assembly: "assembly/",
        PIPELINE_NAME_contig_classification: "classification/assembly/",
        PIPELINE_NAME_read_classification: "classification/reads/",
        PIPELINE_NAME_remapping: "remap/",
        "log_dir": "logs/",
        "OUTD": "output/",
    }

    ################################## MODULES CONTROL

    ACTIONS = {
        "CLEAN": False,
        "QCONTROL": False,
        "ENRICH": True,
        "DEPLETE": False,
        "ASSEMBLE": True,
        "CLASSIFY": True,
        "REMAP": True,
        "PHAGE_DEPL": True,
        "VIRSORT": False,
        "SIFT": True,
    }

    ################################## TECHNOLOGY CONSTANTS

    CONSTANTS_ILLUMINA = {
        "minimum_coverage_threshold": 2,
        "max_output_number": 100,
        "taxid_limit": 100,
        "sift_query": "phage",
        "assembly_contig_min_length": 300,
    }

    CONSTANTS_ONT = {
        "minimum_coverage_threshold": 1,
        "max_output_number": 100,
        "taxid_limit": 100,
        "sift_query": "phage",
        "assembly_contig_min_length": 500,
    }
