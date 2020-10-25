cmake_minimum_required(VERSION 3.17)
project(fa)

set(CMAKE_CXX_STANDARD 20)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pthread -isystem external")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

set(ALIGNER
        src/aligner/Aligner.cpp src/aligner/Aligner.hpp)

set(ALIGNER_METHODS
        src/aligner/methods/hirschberg.cpp src/aligner/methods/hirschberg.hpp)

set(FILE
        src/file/File.hpp src/file/File.cpp src/file/FileManager.hpp src/file/FileManager.cpp)

set(FORMATS
        src/file/formats/FASTA.hpp src/file/formats/SNPS.hpp)

set(GENOME
        src/genome/Genome.cpp src/genome/Genome.hpp src/genome/Reference.cpp src/genome/Reference.hpp
        src/genome/Query.cpp src/genome/Query.hpp)

set(TYPES
        src/types/alignment.hpp)

set(UTILS
        src/utils/searcher.cpp src/utils/searcher.hpp src/utils/threader.hpp)

set(EXCEPTIONS
        src/exceptions/Exceptions.hpp src/exceptions/FileExceptions.hpp src/exceptions/GenomeExceptions.hpp)

add_executable(fa src/main.cpp ${ALIGNER} ${ALIGNER_METHODS} ${FILE} ${FORMATS} ${GENOME} ${TYPES} ${UTILS} ${EXCEPTIONS})