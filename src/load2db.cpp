#include <htslib/sam.h>
#include <bam_sql.h>
#include <bam_utils.h>
#include <iostream>
#include <gperftools/profiler.h>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

const char* convert_to_cstr(const string & s)
{
    return s.c_str();
}

int main(int argc, char** argv) {

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("bam,b", po::value<string>(), "indexed bam")
        ("barcode,c", po::value<string>(), "prefix for barcode file")
        ("db,d", po::value<string>(), "output database")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }


    const char* bam_fname = convert_to_cstr(vm["bam"].as<string>());
    const char* dest_fname = convert_to_cstr(vm["db"].as<string>());
    const char* barcodes = convert_to_cstr(vm["barcode"].as<string>());

    BamDB* bamdb = new BamDB(bam_fname, dest_fname, barcodes);

    create_table(bamdb);
    fill_db(bamdb);

    return 0;
}
