#include <htslib/sam.h>
#include <bam_hash.h>
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
        //("paired-end", po::value<bool>()->default_value(false), "Record matepair information")
        ("barcode", po::value<string>(), "prefix of barcode file")
        ("db", po::value<string>())
        ("anno", po::value<string>())
        ("bam", po::value<string>())
    ;

    po::positional_options_description positionalOptions;
    positionalOptions.add("barcode", 1);
    positionalOptions.add("db", 1);
    positionalOptions.add("anno", 1);
    positionalOptions.add("bam", 1);

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
             << "  db : output database" << endl
             << "  bam : BAM to be parsed" << endl;
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

    const char* bam_fname = convert_to_cstr(vm["bam"].as<string>());
    const char* dest_fname = convert_to_cstr(vm["db"].as<string>());
    const char* anno_fname = convert_to_cstr(vm["anno"].as<string>());
    const char* barcodes = convert_to_cstr(vm["barcode"].as<string>());


    int umi_length = vm["umi-length"].as<int>();
    int bc_min_qual = vm["bc-min-qual"].as<int>();


    BamHash * bamhash = new BamHash(bam_fname, anno_fname, dest_fname, barcodes, umi_length, bc_min_qual);

    hash_annotation(bamhash);
    bamhash->hash_reads();
    bamhash->print_results();
    delete bamhash;
    return 0;
}
