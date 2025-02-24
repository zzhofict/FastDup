/*
Description: Markduplicate需要用到的一些参数

Copyright : All right reserved by ICT

Author : Zhang Zhonghai
Date : 2023/10/23
*/
#pragma once

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace nsmd {
/* How strict to be when reading a SAM or BAM, beyond bare minimum validation.
 */
enum ValidationStringency {
    /**
     * Do the right thing, throw an exception if something looks wrong.
     */
    STRICT,
    /**
     * Emit warnings but keep going if possible.
     */
    LENIENT,
    /**
     * Like LENIENT, only don't emit warning messages.
     */
    SILENT,

    DEFAULT_STRINGENCY = SILENT
};

/**
 * Enum used to control how duplicates are flagged in the DT optional tag on
 * each read.
 */
enum DuplicateTaggingPolicy { DontTag, OpticalOnly, All };

/* 排序的方式 */
enum SortOrder {
    unsorted,
    queryname,
    coordinate,
    duplicate,  // NB: this is not in the SAM spec!
    unknown
};

/* 计算reads分数的方式（比那个read得分更高） */
enum ScoringStrategy { SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM };

/* 索引文件的格式 （bai或者csi） */
enum IndexFormat { BAI, CSI };
}  // namespace nsmd

/* markduplicate 需要的参数*/
struct MarkDupsArg {
    string INPUT_FILE;  // input bam filename

    string OUTPUT_FILE;  // output bam filename

    int NUM_THREADS = 1;

    size_t MAX_MEM = ((size_t)1) << 30; // << 31  // 最小2G

    bool DUPLEX_IO = true; // 同时读写

    /**
     * The optional attribute in SAM/BAM/CRAM files used to store the duplicate type.
     */
    string DUPLICATE_TYPE_TAG = "DT";
    /**
     * The duplicate type tag value for duplicate type: library.
     */
    string DUPLICATE_TYPE_LIBRARY = "LB";
    /**
     * The duplicate type tag value for duplicate type: sequencing (optical & pad-hopping, or "co-localized").
     */
    string DUPLICATE_TYPE_SEQUENCING = "SQ";
    /**
     * The attribute in the SAM/BAM file used to store which read was selected as representative out of a duplicate set
     */
    string DUPLICATE_SET_INDEX_TAG = "DI";
    /**
     * The attribute in the SAM/BAM file used to store the size of a duplicate set
     */
    string DUPLICATE_SET_SIZE_TAG = "DS";

    /* OpticalDuplicateFinder */
    int DEFAULT_OPTICAL_DUPLICATE_DISTANCE = 100;
    int DEFAULT_BIG_DUPLICATE_SET_SIZE = 1000;
    int DEFAULT_MAX_DUPLICATE_SET_SIZE =
        300000;  // larger than this number will generate over 100 billion comparisons in the n^2 algorithm below

    /**
     * If more than this many sequences in SAM file, don't spill to disk because there will not
     * be enough file handles.
     */

    /* "This option is obsolete. ReadEnds will always be spilled to disk." */
    int MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP = 50000;

    /* "Maximum number of file handles to keep open when spilling read ends to disk. " +
                    "Set this number a little lower than the per-process maximum number of file that may be open. " +
                    "This number can be found by executing the 'ulimit -n' command on a Unix system." */
    int MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 8000;

    /* "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by " +
                    "some of the sorting collections.  If you are running out of memory, try reducing this number." */
    double SORTING_COLLECTION_SIZE_RATIO = 0.25;

    /* "Barcode SAM tag (ex. BC for 10X Genomics)", optional = true */
    string BARCODE_TAG = "";

    /* "Read one barcode SAM tag (ex. BX for 10X Genomics)", optional = true */
    string READ_ONE_BARCODE_TAG = "";

    /* "Read two barcode SAM tag (ex. BX for 10X Genomics)", optional = true */
    string READ_TWO_BARCODE_TAG = "";

    /* "If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), " +
                     "indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two
       " + "reads map to the same portion of the reference only one of which is marked as duplicate. The second " +
                     "tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which
       the " + "record belongs. This identifier is the index-in-file of the representative read that was selected out "
       + "of the duplicate set.", optional = true) */
    bool TAG_DUPLICATE_SET_MEMBERS = false;

    /* "If true remove 'optical' duplicates and other duplicates that appear to have arisen from the " +
                    "sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false.
       " + "If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored.") */
    bool REMOVE_SEQUENCING_DUPLICATES = false;

    /* "Determines how duplicate types are recorded in the DT optional attribute.") */
    nsmd::DuplicateTaggingPolicy TAGGING_POLICY = nsmd::DuplicateTaggingPolicy::DontTag;

    /* "Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this tag.  Default true")
     */
    bool CLEAR_DT = true;

    /* "Treat UMIs as being duplex stranded.  This option requires that the UMI consist of two equal length " +
                    "strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered duplicates if, in
       addition to standard " + "definition, have identical normalized UMIs.  A UMI from the 'bottom' strand is
       normalized by swapping its content " + "around the hyphen (eg. ATC-GTC becomes GTC-ATC).  A UMI from the 'top'
       strand is already normalized as it is. " + "Both reads from a read pair considered top strand if the read 1
       unclipped 5' coordinate is less than the read " + "2 unclipped 5' coordinate. All chimeric reads and read
       fragments are treated as having come from the top strand. " + "With this option is it required that the
       BARCODE_TAG hold non-normalized UMIs. Default false.") */
    bool DUPLEX_UMI = false;

    /* "SAM tag to uniquely identify the molecule from which a read was derived.  Use of this option requires that " +
                    "the BARCODE_TAG option be set to a non null value.  Default null.",
              optional = true) */
    string MOLECULAR_IDENTIFIER_TAG = "";

    /* 继承自 AbstractMarkDuplicatesCommandLineProgram 的参数*/
    /* "File to write duplication metrics to" */
    string METRICS_FILE;

    /* "If true do not write duplicates to the output file instead of writing them with appropriate flags set." */
    bool REMOVE_DUPLICATES = false;

    /* "If true, assume that the input file is coordinate sorted even if the header says otherwise. " +
                    "Deprecated, used ASSUME_SORT_ORDER=coordinate instead." mutex = {"ASSUME_SORT_ORDER"} */
    bool ASSUME_SORTED = false;

    /* "If not null, assume that the input file has this order even if the header says otherwise.",
              optional = true, mutex = {"ASSUME_SORTED"} */
    nsmd::SortOrder ASSUME_SORT_ORDER = nsmd::SortOrder::unsorted;

    /* "The scoring strategy for choosing the non-duplicate among candidates." */
    nsmd::ScoringStrategy DUPLICATE_SCORING_STRATEGY = nsmd::ScoringStrategy::SUM_OF_BASE_QUALITIES;
    // nsmd::ScoringStrategy DUPLICATE_SCORING_STRATEGY = nsmd::ScoringStrategy::TOTAL_MAPPED_REFERENCE_LENGTH;

    /* "The program record ID for the @PG record(s) created by this program. Set to null to disable " +
                    "PG record creation.  This string may have a suffix appended to avoid collision with other " +
                    "program record IDs.",
              optional = true */
    string PROGRAM_RECORD_ID = "FastDup";

    /* "Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.",
              optional = true */
    string PROGRAM_GROUP_VERSION;

    /* "Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.",
              optional = true */
    string PROGRAM_GROUP_COMMAND_LINE;

    /* "Value of PN tag of PG record to be created." */
    string PROGRAM_GROUP_NAME = "FastDup";

    /* "Comment(s) to include in the output file's header.",
              optional = true */
    vector<string> COMMENT;

    /* 继承自 AbstractOpticalDuplicateFinderCommandLineProgram 的参数 */

    /* "MarkDuplicates can use the tile and cluster positions to estimate the rate of optical duplication " +
            "in addition to the dominant source of duplication, PCR, to provide a more accurate estimation of library
       size. " + "By default (with no READ_NAME_REGEX specified), MarkDuplicates will attempt to extract coordinates " +
            "using a split on ':' (see Note below).  " +
            "Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. " +
            "Note that without optical duplicate counts, library size estimation will be less accurate. " +
            "If the read name does not follow a standard Illumina colon-separation convention, but does contain tile and
       x,y coordinates, " + "a regular expression can be specified to extract three variables: tile/region, x coordinate
       and y coordinate from a read name. " + "The regular expression must contain three capture groups for the three
       variables, in order. " + "It must match the entire read name. " + "  e.g. if field names were separated by
       semi-colon (';') this example regex could be specified " + " (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ "
       + "Note that if no READ_NAME_REGEX is specified, the read name is split on ':'. " + "  For 5 element names, the
       3rd, 4th and 5th elements are assumed to be tile, x and y values. " + "  For 7 element names (CASAVA 1.8), the
       5th, 6th, and 7th elements are assumed to be tile, x and y values.", optional = true */
    string READ_NAME_REGEX = "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$";
    bool CHECK_OPTICAL_DUP = true;

    /* "The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default " +
                    "is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell
       models, 2500 is more" + "appropriate. For other platforms and models, users should experiment to find what works
       best." */
    int OPTICAL_DUPLICATE_PIXEL_DISTANCE = DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

    /* "This number is the maximum size of a set of duplicate reads for which we will attempt to determine " +
                    "which are optical duplicates.  Please be aware that if you raise this value too high and do
       encounter a very " + "large set of duplicate reads, it will severely affect the runtime of this tool.  To
       completely disable this check, " + "set the value to -1." */
    long MAX_OPTICAL_DUPLICATE_SET_SIZE = DEFAULT_MAX_DUPLICATE_SET_SIZE;

    /* 继承自 CommandLineProgram 的参数*/

    /* "Whether to suppress job-summary info on System.err.", common = true */
    bool QUIET = false;

    /* "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT " +
            "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) "
       + "do not otherwise need to be decoded.", common=true */
    nsmd::ValidationStringency VALIDATION_STRINGENCY = nsmd::ValidationStringency::DEFAULT_STRINGENCY;

    /* "Compression level for all compressed files created (e.g. BAM and VCF).", common = true */
    int COMPRESSION_LEVEL = 5;

    /* "When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling
       to disk. " + "Increasing this number reduces the number of file handles needed to sort the file, and increases
       the amount of RAM needed.", optional = true, common = true */
    int MAX_RECORDS_IN_RAM = 500000;

    /* "Whether to create an index when writing VCF or coordinate sorted BAM output.", common = true */
    bool CREATE_INDEX = false;

    nsmd::IndexFormat INDEX_FORMAT = nsmd::IndexFormat::BAI;

    /* "Whether to create an MD5 digest for any BAM or FASTQ files created.  ", common = true */
    bool CREATE_MD5_FILE = false;

    /* Add PG tag to each read in a SAM or BAM (PGTagArgumentCollection)*/
    bool ADD_PG_TAG_TO_READS = true;

    // 命令行字符串
    string CLI_STR;

    // 开始运行时间
    string START_TIME;
};