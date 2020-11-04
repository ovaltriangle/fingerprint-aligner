#include "aligner/Aligner.hpp"
#include "exceptions/Exceptions.hpp"

#include <iostream>
#include <memory>

#include <unistd.h>


#include "../external/cxxopts/cxxopts.hpp"

int main(const int argc, const char** argv) {
    const std::string version {"1.0"};

    const std::string default_out_seq {"out.fa"};
    const std::string default_out_snps {"out.snps"};
    const std::string default_threads {"1"};
    const std::string default_p_thresh {"60"};
    const std::string default_s_thresh {"0.9"};
    const std::string default_weights {"-2,-2,-1,1"};

    cxxopts::Options options("fa", "Aligns genomes using a reference fingerprint.");

    options.add_options()
            ("r,reference", "The file containing the reference genomes to be used.", cxxopts::value<std::string>())
            ("q,query", "The file containing the query genomes(s).", cxxopts::value<std::string>())
            ("o,output", "Sequence file output. If not specified, '" + default_out_seq + "' will be used.",
             cxxopts::value<std::string>()->default_value(default_out_seq))
//            ("n,snp", "SNPs file output. If not specified, '" + default_out_snps + "' will be used.",
//            cxxopts::value<std::string>()->default_value(default_out_snps))
            ("t,threads", "The number of threads to use (>1 for multithreading).",
             cxxopts::value<int>()->default_value(default_threads))
            ("p,protein-threshold", "The threshold for the number of nucleotides before a protein is considered so.",
             cxxopts::value<std::size_t>()->default_value(default_p_thresh))
            ("s,similarity-threshold", "Minimum percentage of identity when comparing two sequences.",
             cxxopts::value<double>()->default_value(default_s_thresh))
            ("w,weights", "The weights to be used in the Hirschberg's NW score calculation"
                          "(insertion, deletion, substitution and match scores). Separate the values with a comma.",
            cxxopts::value< std::vector<long> >()->default_value(default_weights))
            // ("m,method", "The method number to use.", cxxopts::value<unsigned short>()->default_value("0"))
            ("h,help", "Print this message and exit.")
            ;

    try {
        auto result {options.parse(argc, argv)};

        if (result.count("help") or result.arguments().empty()) {
            std::cout << options.help() << '\n';
            exit(0);
        }

        if (result.count("r") != 1) {
            std::cerr << "Be sure to specify one reference files.\n\n";
            std::cerr << options.help() << std::endl;
            exit(1);
        }

        if (result.count("q") != 1) {
            std::cerr << "Be sure to specify one query files.\n\n";
            std::cerr << options.help() << std::endl;
            exit(1);
        }

        int threads_num {0};
        if (result.count("t")) {
            threads_num = result["t"].as<int>();
            threads_num = (threads_num > 1 ? threads_num - 1 : 0);
        }

        auto weights_pre {result["w"].as<std::vector<long>>()};
        WeightTable weights {weights_pre.at(0), weights_pre.at(1),
                             weights_pre.at(2), weights_pre.at(3)};

        std::string reference {result["r"].as<std::string>()};
        std::string query {result["q"].as<std::string>()};
        std::string output {result["o"].as<std::string>()};
        //std::string snps {result["n"].as<std::string>()};
        std::string snps {"mock.fasta"};
        // threads_num already has the thread number
        std::size_t p_threshold {result["p"].as<std::size_t>()};
        double s_threshold {result["s"].as<double>()};
        // weights have also been preprocessed
        // const unsigned short *method {&result["m"].as<unsigned short >()};

        if (reference == query) {
            throw RQSameFile();
        }

        // checkers for overwriting files

        try {
            /// ALIGNER
            std::cout << "Starting Fingerprint Aligner v" + version + '\n';
            auto aligner {std::make_unique<Aligner>(reference, query, output, snps, p_threshold, s_threshold, weights)};
            aligner->perform_alignment(threads_num);
            std::cout << "Success! Everything went according to plan.\n";
        } catch (std::exception &e) {
            std::cerr << "An exception occurred while handling your request:\n\t" << e.what();
        }

    } catch (cxxopts::OptionParseException &e) {
        std::cerr << "There was an error while parsing the arguments:\n\t" << e.what() << '\n';
    } catch (std::exception &e) {
        std::cerr << "An exception occurred while parsing your arguments:\n\t" << e.what();
    }
}
