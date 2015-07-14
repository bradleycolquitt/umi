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

const char* convert_to_cstr(const string & s)
{
    return s.c_str();
}

int main(int argc, char** argv) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("umi-length,u", po::value<int>()->default_value(8), "length of UMI sequence")
        ("bc-min-qual,q", po::value<int>()->default_value(20), "minimum quality score for barcode")
        ("i5", po::value<int>()->default_value(0), "i5 index sequence")
        ("i7", po::value<int>()->default_value(0), "i7 index sequence")
        ("print-to-txt,p", po::bool_switch()->default_value(false), "export results to plain text instead of sqlite3 db")
        ("barcode", po::value<string>(), "prefix of barcode file")
        ("anno", po::value<string>())
      //  ("bam", po::value<string>())
        ("fastq", po::value<string>())
        ("outfile", po::value<string>())
    ;

    po::positional_options_description positionalOptions;
    positionalOptions.add("barcode", 1);
    positionalOptions.add("anno", 1);
    //positionalOptions.add("bam", 1);
    positionalOptions.add("fastq", 1);
    positionalOptions.add("outfile", 1);

    po::variables_map vm;

    try
    {
      po::store(po::command_line_parser(argc, argv).options(desc)
                  .positional(positionalOptions).run(),
                vm); // throws on error

      /** --help option
       */
      if ( vm.count("help")  )
      {
        std::cout << endl
                  << "Parses and loads reads from a BAM file to an SQLite3 DB." << endl
                  << "Expects following read structure: "
                  << endl << endl
                  << "\t[UMI][Barcode][Variable number of G][Sequence]"
                  << endl << endl;
        cout << "Usage: load2db [options] barcode_prefix db bam" << endl << endl
             << "Required options:" << endl
             << "  barcode_prefix : prefix of file found in barcode directory" << endl
             << "  anno : .featureCounts file"  << endl
             //<< "  bam : BAM to be parsed" << endl
             << "  outfile : output file" << endl;
        cout << desc << endl;
        return SUCCESS;
      }

      po::notify(vm); // throws on error, so do after help in case
                      // there are any problems
    }
    catch(boost::program_options::required_option& e)
    {
      std::cerr << "ERROR: " << e.what() << ", required." << std::endl << std::endl;
      return ERROR_IN_COMMAND_LINE;
    }
    catch(boost::program_options::error& e)
    {
      std::cerr << "ERROR: " << e.what() << ", optional." << std::endl << std::endl;
      return ERROR_IN_COMMAND_LINE;
    }

    //const char* bam_fname = convert_to_cstr(vm["bam"].as<string>());
    const char* fastq_fname = convert_to_cstr(vm["fastq"].as<string>());
    const char* dest_fname = convert_to_cstr(vm["outfile"].as<string>());
    const char* anno_fname = convert_to_cstr(vm["anno"].as<string>());
    const char* barcodes = convert_to_cstr(vm["barcode"].as<string>());


    int umi_length = vm["umi-length"].as<int>();
    int bc_min_qual = vm["bc-min-qual"].as<int>();
    int i5 = vm["i5"].as<int>();
    int i7 = vm["i7"].as<int>();
    bool to_txt = vm["print-to-txt"].as<bool>();

    //BamHash * bamhash = new BamHash(bam_fname, fastq_fname, anno_fname, dest_fname, barcodes, umi_length, bc_min_qual, i5, i7, to_txt);
    BamHash * bamhash = new BamHash(fastq_fname, anno_fname, dest_fname, barcodes, umi_length, bc_min_qual, i5, i7, to_txt);

    hash_annotation(bamhash);
    bamhash->hash_reads();

    cout << "----Outputting results----" << endl;
    if (to_txt == 1)
    {
        bamhash->print_results();
    }
    else
    {
        bamhash->write_to_db();
    }

    delete bamhash;
    return 0;
}
