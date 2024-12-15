#include <spdlog/cfg/env.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <argparse/argparse.hpp>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "markdup/md_args.h"
#include "util/profiling.h"
#include "fastdup_version.h"

namespace nsgv {
extern MarkDupsArg gMdArg;
};

int MarkDuplicates();

int main(int argc, char *argv[]) {
    // init log
    spdlog::set_default_logger(spdlog::stderr_color_st("fastdup"));
    spdlog::cfg::load_env_levels();
    
    // init arg parser
    argparse::ArgumentParser program(nsgv::gMdArg.PROGRAM_RECORD_ID, FASTDUP_VERSION, argparse::default_arguments::none);
    program.add_description(
        "Identifies duplicate reads. This tool locates and tags duplicate reads in a coordinate ordered SAM or BAM "
        "file.\nUse the same algorithm as picard MarkDuplicates and output identical results.\n"
        "Use spdlog as log tool and the default level is 'info'.");

    program.add_argument("--input")
        .help("Input file. One coordinate ordered SAM or BAM file.")
        .metavar("<INPUT>")
        .required();

    program.add_argument("--metrics")
        .help("Metrics file. File to write duplication metrics to.")
        .metavar("<METRICS>")
        .required();

    program.add_argument("--output")
        .help("Output file. SAM or BAM file to write marked records to.")
        .metavar("<OUTPUT>")
        .required();

    program.add_argument("--num-threads")
        .help("Number of threads to use.")
        .scan<'i', int>()
        .default_value(1)
        .nargs(1)
        .metavar("<NUM_THREADS>");

    program.add_argument("--none-duplex-io")
        .help("Do not use writing-while-reading mode.")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--create-index")
        .help("Whether to create an index when writing coordinate sorted BAM output.")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--index-format")
        .help("Format for bam index file. Possible values: {BAI, CSI}")
        .default_value(std::string("BAI"))
        .choices("BAI", "CSI")
        .nargs(1)
        .metavar("<IndexFormat>");

    program.add_argument("--remove-duplicates")
        .help("If true do not write duplicates to the output file instead of writing them with appropriate flags set.")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--duplicate-scoring-strategy")
        .help(
            "The scoring strategy for choosing the non-duplicate among candidates. "
            "Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}")
        .default_value(std::string("SUM_OF_BASE_QUALITIES"))
        .choices("SUM_OF_BASE_QUALITIES", "TOTAL_MAPPED_REFERENCE_LENGTH", "RANDOM")
        .nargs(1)
        .metavar("<ScoringStrategy>");

    program.add_argument("--optical-duplicate-pixel-distance")
        .help(
            "\nThe maximum offset between two duplicate clusters in order to consider them optical"
            "duplicates. The default is appropriate for unpatterned versions of the Illumina platform."
            "For the patterned flowcell models, 2500 is moreappropriate. For other platforms and"
            "models, users should experiment to find what works best.")
        .scan<'i', int>()
        .default_value(100)
        .nargs(1)
        .metavar("<Integer>");

    program.add_argument("--read-name-regex")
        .help(
            "\nMarkDuplicates can use the tile and cluster positions to estimate the rate of optical "
            "duplication in addition to the dominant source of duplication, PCR, to provide a more "
            "accurate estimation of library size. By default (with no READ_NAME_REGEX specified), "
            "MarkDuplicates will attempt to extract coordinates using a split on ':' (see Note below). "
            "Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. Note that without "
            "optical duplicate counts, library size estimation will be less accurate. If the read name "
            "does not follow a standard Illumina colon-separation convention, but does contain tile and "
            "x,y coordinates, a regular expression can be specified to extract three variables: "
            "tile/region, x coordinate and y coordinate from a read name. The regular expression must "
            "contain three capture groups for the three variables, in order. It must match the entire "
            "read name.   e.g. if field names were separated by semi-colon (';') this example regex "
            "could be specified      (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ Note that if no "
            "READ_NAME_REGEX is specified, the read name is split on ':'.   For 5 element names, the "
            "3rd, 4th and 5th elements are assumed to be tile, x and y values.   For 7 element names "
            "(CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. "
            "Default value: <optimized capture of last three ':' separated fields as numeric values>.")
        .default_value(std::string("(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$"))
        .nargs(1)
        .metavar("<ReadNameRegex>");

    program.add_argument("--tag-duplicate-set-members")
        .help(
            "\nIf a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG"
            "(DS), indicates the size of the duplicate set. The smallest possible DS value is 2 which"
            "occurs when two reads map to the same portion of the reference only one of which is marked"
            "as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier"
            "for the duplicate set to which the record belongs. This identifier is the index-in-file of"
            "the representative read that was selected out of the duplicate set.")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--tagging-policy")
        .help(
            "Determines how duplicate types are recorded in the DT optional attribute. Possible values: {DontTag, "
            "OpticalOnly, All}.")
        .default_value(std::string("DontTag"))
        .choices("DontTag", "OpticalOnly", "All")
        .nargs(1)
        .metavar("<DuplicateTaggingPolicy>");

    // add help and version args
    program.add_argument("-h", "--help")
        .action([&](const auto & /*unused*/) {
            std::cout << program.help().str();
            std::exit(0);
        })
        .default_value(false)
        .help("shows help message and exits")
        .implicit_value(true)
        .nargs(0);

    program.add_argument("-v", "--version")
        .action([&](const auto & /*unused*/) {
            std::cout << FASTDUP_VERSION << std::endl;
            std::exit(0);
        })
        .default_value(false)
        .help("prints version information and exits")
        .implicit_value(true)
        .nargs(0);

    // std::cout << program << std::endl;

    try {
        program.parse_args(argc, argv);
        nsgv::gMdArg.INPUT_FILE = program.get("--input");
        nsgv::gMdArg.OUTPUT_FILE = program.get("--output");
        nsgv::gMdArg.METRICS_FILE = program.get("--metrics");
        nsgv::gMdArg.NUM_THREADS = program.get<int>("--num-threads");
        nsgv::gMdArg.DUPLEX_IO = !program.get<bool>("--none-duplex-io");
        nsgv::gMdArg.CREATE_INDEX = program.get<bool>("--create-index");
        
        nsgv::gMdArg.INDEX_FORMAT =
            program.get("--index-format") == "BAI" ? nsmd::IndexFormat::BAI : nsmd::IndexFormat::CSI;
        nsgv::gMdArg.REMOVE_DUPLICATES = program.get<bool>("--remove-duplicates");
        
        std::map<std::string, nsmd::ScoringStrategy> scoring_strategy_args = {
            {"SUM_OF_BASE_QUALITIES", nsmd::ScoringStrategy::SUM_OF_BASE_QUALITIES},
            {"TOTAL_MAPPED_REFERENCE_LENGTH", nsmd::ScoringStrategy::TOTAL_MAPPED_REFERENCE_LENGTH},
            {"RANDOM", nsmd::ScoringStrategy::RANDOM}};
        nsgv::gMdArg.DUPLICATE_SCORING_STRATEGY = scoring_strategy_args[program.get("--duplicate-scoring-strategy")];
        
        nsgv::gMdArg.OPTICAL_DUPLICATE_PIXEL_DISTANCE = program.get<int>("--optical-duplicate-pixel-distance");
        
        std::set<string> all_nulls = {"null", "Null", "NULL"};
        nsgv::gMdArg.READ_NAME_REGEX =
            all_nulls.find(program.get("--read-name-regex")) != all_nulls.end() ? "" : program.get("--read-name-regex");
        if (nsgv::gMdArg.READ_NAME_REGEX.empty())
            nsgv::gMdArg.CHECK_OPTICAL_DUP = false;
        
        nsgv::gMdArg.TAG_DUPLICATE_SET_MEMBERS = program.get<bool>("--tag-duplicate-set-members");
        
        std::map<std::string, nsmd::DuplicateTaggingPolicy> tagging_policy_args = {
            {"DontTag", nsmd::DuplicateTaggingPolicy::DontTag},
            {"OpticalOnly", nsmd::DuplicateTaggingPolicy::OpticalOnly},
            {"All", nsmd::DuplicateTaggingPolicy::All}};
        nsgv::gMdArg.TAGGING_POLICY = tagging_policy_args[program.get("--tagging-policy")];

    } catch (const std::exception &err) {
        spdlog::error(err.what());
        return 1;
    }

    spdlog::info("fast markduplicates start");
    MarkDuplicates();
    spdlog::info("fast markduplicates end");

    DisplayProfiling(1);

    return 0;
}