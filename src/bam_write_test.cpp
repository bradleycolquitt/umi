#include <htslib/sam.h>
#include <bam_hash2.h>
#include <bamrecord.h>
#include <bam_utils.h>
#include <iostream>
#include <gperftools/profiler.h>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

namespace
{
const size_t ERROR_IN_COMMAND_LINE = 1;
const size_t SUCCESS = 0;
const size_t ERROR_UNHANDLED_EXCEPTION = 2;
} // namespace

int main(int argc, char** argv) {
    const char* bam_fname = "../test_files/test.bam";
    string out_fname = "../test_files/test_out.bam";
    const char* fastq_fname = "../test_files/test.fastq.gz";
    const char* anno_fname = "../test_files/test.featureCounts";
    const char* barcodes = "bc24";
    const char* dest_fname = "../test_files/test.db";
    //const char* bam_fname = convert_to_cstr(vm["bam"].as<string>());
    // const char* fastq_fname = convert_to_cstr(vm["fastq"].as<string>());
    // const char* dest_fname = convert_to_cstr(vm["outfile"].as<string>());
    // const char* anno_fname = convert_to_cstr(vm["anno"].as<string>());
    // const char* barcodes = convert_to_cstr(vm["barcode"].as<string>());

    int umi_length = 25;
    int bc_min_qual = 10;
    int i5 = 0;
    int i7 = 1;
    bool to_txt = false;
    BamHash * bamhash = new BamHash(bam_fname, fastq_fname, anno_fname, dest_fname, barcodes, umi_length, bc_min_qual, i5, i7, to_txt);
    //BamHash * bamhash = new BamHash(fastq_fname, anno_fname, dest_fname, barcodes, umi_length, bc_min_qual, i5, i7, to_txt);

    //bamhash->write_bam(out_fname);
    //hash_bam(bamhash);
    bamhash->split_bam();
    //hash_annotation(bamhash);
    //bamhash->hash_reads();

    delete bamhash;
    return 0;
}
