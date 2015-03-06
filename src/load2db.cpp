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
        ("help,h", "produce help message")
        ("umi_positions,u", po::value<vector<int> >()->multitoken(), "start and end positions of UMI sequence")
        ("bam", po::value<string>()->required(), "indexed bam")
        ("barcode", po::value<string>()->required(), "prefix for barcode file")
        ("db", po::value<string>()->required(), "output database")
    ;

    po::positional_options_description positionalOptions;
    positionalOptions.add("bam", 1);
    positionalOptions.add("barcode", 1);
    positionalOptions.add("db", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc)
                .positional(positionalOptions).run(),
                          vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        //return 1;
    }

    po::notify(vm);


    const char* bam_fname = convert_to_cstr(vm["bam"].as<string>());
    const char* dest_fname = convert_to_cstr(vm["db"].as<string>());
    const char* barcodes = convert_to_cstr(vm["barcode"].as<string>());

    vector<int> umi_pos;
    if ( ! vm.count("umi_positions") ) {
        umi_pos.push_back(0);
        umi_pos.push_back(4);
    } else {
        umi_pos = vm["umi_positions"].as<vector<int> >();
    }
    BamDB* bamdb = new BamDB(bam_fname, dest_fname, barcodes, umi_pos[0], umi_pos[1]);

    create_table(bamdb);
    fill_db(bamdb);
    bamdb->create_reftable();
    create_index(bamdb);

    return 0;
}
