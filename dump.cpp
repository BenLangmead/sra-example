//
// Created by Ben Langmead on 10/3/18.
//

#include <iostream>
#include <string>
#include <utility>
#include <stdexcept>
#include "dump.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::pair;

#if defined(_TTHREAD_WIN32_)
#define SLEEP(x) Sleep(x)
#else
#define SLEEP(x) do { \
	const static timespec ts_tmp_ = {0, 1000000 * x}; \
	nanosleep(&ts_tmp_, NULL); \
} while(false)
#endif

static const uint64_t buffer_size_per_thread = 4096;

struct SRA_Data {
    uint64_t read_pos;
    uint64_t write_pos;
    uint64_t buffer_size;
    bool     done;

    ngs::ReadIterator* sra_it;

    SRA_Data() {
        read_pos = 0;
        write_pos = 0;
        buffer_size = buffer_size_per_thread;
        done = false;
        sra_it = NULL;
    }

    bool isFull() {
        assert(read_pos <= write_pos);
        assert(read_pos + buffer_size >= write_pos);
        return read_pos + buffer_size <= write_pos;
    }
};

void downloader(SRA_Data* sra_data) {
    assert(sra_data != NULL);
    ngs::ReadIterator* sra_it = sra_data->sra_it;
    assert(sra_it != NULL);

    while(!sra_data->done) {
        while(sra_data->isFull()) { SLEEP(1); }
        bool exception_thrown = false;
        try {
            if(!sra_it->nextRead() || !sra_it->nextFragment()) {
                return;
            }
            // Read the name out of the first field
            ngs::StringRef rname = sra_it->getReadId();
            ngs::StringRef ra_seq = sra_it->getFragmentBases();
            ngs::StringRef ra_qual = sra_it->getFragmentQualities();
            if(!sra_it->nextFragment()) {
                // done
                //cout << "Read: " << rname << ":" << ra_seq << ":" << ra_qual << endl;
            } else {
                // rb.name = ra.name;
                ngs::StringRef rb_seq = sra_it->getFragmentBases();
                ngs::StringRef rb_qual = sra_it->getFragmentQualities();
                //cout << "Read: " << rname << ":" << ra_seq << ":" << ra_qual << ":" << rb_seq << ":" << rb_qual << endl;
            }
            sra_data->write_pos++;
        } catch(ngs::ErrorMsg & x) {
            cerr << x.toString () << endl;
            exception_thrown = true;
        } catch(std::exception & x) {
            cerr << x.what () << endl;
            exception_thrown = true;
        } catch(...) {
            cerr << "unknown exception\n";
            exception_thrown = true;
        }

        if(exception_thrown) {
            sra_data->done = true;
            cerr << "An error happened while fetching SRA reads. Please rerun HISAT2. You may want to disable the SRA cache if you didn't (see the instructions at https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration).\n";
            exit(1);
        }
    }
}

void open(const string& sra_acc) {
    string version = "sra-example-btl";
    ncbi::NGS::setAppVersionString(version);
    assert(!sra_acc.empty());
    try {
        // open requested accession using SRA implementation of the API
        ngs::ReadCollection sra_run = ncbi::NGS::openReadCollection(sra_acc);

        // compute window to iterate through
        size_t MAX_ROW = sra_run.getReadCount();
        cerr << "MAX_ROW = " << MAX_ROW << endl;
        ngs::ReadIterator* sra_it;
        sra_it = new ngs::ReadIterator(sra_run.getReadRange(1, MAX_ROW, ngs::Read::all));

        // create a buffer for SRA data
        SRA_Data *sra_data = new SRA_Data;
        sra_data->sra_it = sra_it;
        sra_data->buffer_size = 1 * buffer_size_per_thread;

        // ready to rock!
    } catch(...) {
        cerr << "Warning: Could not access \"" << sra_acc << "\" for reading; skipping..." << endl;
    }
}

int main(int argc, char **argv) {
    if(argc < 2) {
        cerr << "Must specify accession as argument" << endl;
    }
    cerr << "Processing accession \"" << argv[1] << "\"" << endl;
    string acc(argv[1]);
    open(acc);
    return 0;
}

