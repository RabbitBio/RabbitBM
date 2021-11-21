#include "barcodeToPositionMulti.h"

#include "../lib/gzip_decompress.hpp" //FIXME

#include "prog_util.h"
#include <fstream>

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <unistd.h>
#include <utime.h>
#include "pigz.h"


double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}


double
PUGZGetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

struct options {
    bool count_lines;
    unsigned nthreads;
};

static const tchar *const optstring = T(":hnlt:V");

static void
show_usage(FILE *fp) {
    fprintf(fp,
            "Usage: %" TS " [-l] [-t n] FILE...\n"
            "Decompress the specified FILEs.\n"
            "\n"
            "Options:\n"
            "  -l        count line instead of content to standard output\n"
            "  -t n      use n threads\n"
            "  -h        print this help\n"
            "  -V        show version and legal information\n",
            program_invocation_name);
}

static void
show_version(void) {
    printf("Pugz parallel gzip decompression program\n"
           "Copyright 2019 Maël Kerbiriou, Rayan Chikhi\n"
           "\n"
           "Based on:\n"
           "gzip compression program v"
    LIBDEFLATE_VERSION_STRING
    "\n"
    "Copyright 2016 Eric Biggers\n"
    "\n"
    "This program is free software which may be modified and/or redistributed\n"
    "under the terms of the MIT license.  There is NO WARRANTY, to the extent\n"
    "permitted by law.  See the COPYING file for details.\n");
}

static int
stat_file(struct file_stream *in, stat_t *stbuf, bool allow_hard_links) {
    if (tfstat(in->fd, stbuf) != 0) {
        msg("%" TS ": unable to stat file", in->name);
        return -1;
    }

    if (!S_ISREG(stbuf->st_mode) && !in->is_standard_stream) {
        msg("%" TS " is %s -- skipping", in->name, S_ISDIR(stbuf->st_mode) ? "a directory" : "not a regular file");
        return -2;
    }

    if (stbuf->st_nlink > 1 && !allow_hard_links) {
        msg("%" TS " has multiple hard links -- skipping "
            "(use -f to process anyway)",
            in->name);
        return -2;
    }

    return 0;
}


void BarcodeToPositionMulti::getMbpmap() {
    mbpmap = new BarcodePositionMap(mOptions);

}


BarcodeToPositionMulti::BarcodeToPositionMulti(Options *opt) {
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mOutStream = NULL;
    mZipFile = NULL;
    mWriter = NULL;
    mUnmappedWriter = NULL;
    bool isSeq500 = opt->isSeq500;
//    mbpmap = new BarcodePositionMap(opt);
//    printf("test4 val is %d\n", mbpmap->GetHashHead()[109547259]);

    //barcodeProcessor = new BarcodeProcessor(opt, &mbpmap->bpmap);
    if (!mOptions->transBarcodeToPos.fixedSequence.empty() || !mOptions->transBarcodeToPos.fixedSequenceFile.empty()) {
        fixedFilter = new FixedFilter(opt);
        filterFixedSequence = true;
    }
    if (mOptions->usePugz) {
        pugzQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(256);
        pugzQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(256);
    }
    if (mOptions->usePigz) {
        pigzQueue = new moodycamel::ReaderWriterQueue<std::pair<int, std::pair<char *, int>>>(256);
        pigzLast.first = new char[1 << 24];
        pigzLast.second = 0;
    }
    pugz1Done = 0;
    pugz2Done = 0;
    producerDone = 0;
    writerDone = 0;
    mergeDone = 0;
    if (mOptions->numPro == 1)mergeDone = 1;
//    cout << "mergeDone " << mergeDone << endl;
}

BarcodeToPositionMulti::~BarcodeToPositionMulti() {
    //if (fixedFilter) {
    //	delete fixedFilter;
    //}
    //unordered_map<uint64, Position*>().swap(misBarcodeMap);
}


void BarcodeToPositionMulti::mergeWrite() {
    double t0 = GetTime();
    int endTag = 0;
    char *tmpData = NULL;
    int tmpSize = 0;
    long long totSize = 0;
    int cntEnd = 0;
    while (!endTag) {
        MPI_Recv(&tmpSize, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        MPI_Barrier(MPI_COMM_WORLD);
        if (tmpSize == -1) {
            cntEnd++;
#ifdef PRINT_INFO

            printf("processor 0 get -1 %d %d\n", mOptions->numPro - 1, cntEnd);
#endif
            if (cntEnd == mOptions->numPro - 1)
                endTag = 1;
        } else {
            tmpData = new char[tmpSize + 1];
//            memset(tmpData, 0, (tmpSize + 1) * sizeof(char));

            MPI_Recv(tmpData, tmpSize, MPI_CHAR, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            MPI_Barrier(MPI_COMM_WORLD);
            totSize += tmpSize;
            mOutputMtx.lock();
            mWriter->inputFromMerge(tmpData, tmpSize);
            mOutputMtx.unlock();
        }
    }
#ifdef PRINT_INFO

    printf("processor merge get %lld data\n", totSize);
    printf("processor 0 get data done\n");
    printf("merge done, cost %.4f\n", GetTime() - t0);
#endif
}

void BarcodeToPositionMulti::pigzWrite() {
/*
 argc 9
argv ./pigz
argv -p
argv 16
argv -k
argv -4
argv -f
argv -b
argv -4096
argv p.fq
 */
    int cnt = 9;
    string outFile = mOptions->transBarcodeToPos.out1;

    char **infos = new char *[9];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = mOptions->pigzThread;
//    printf("th num is %d\n", th_num);
    string th_num_s = to_string(th_num);
//    printf("th num s is %s\n", th_num_s.c_str());
//    printf("th num s len is %d\n", th_num_s.length());

    infos[2] = new char[th_num_s.length() + 1];
    memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";
    infos[4] = "-4";
    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    string out_file = mOptions->transBarcodeToPos.out1;
//    printf("th out_file is %s\n", out_file.c_str());
//    printf("th out_file len is %d\n", out_file.length());
    infos[8] = new char[out_file.length() + 1];
    memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';

    main_pigz(cnt, infos, pigzQueue, &writerDone, pigzLast);
}

bool BarcodeToPositionMulti::process() {
    auto t0 = GetTime();

    initOutput();
    initPackRepositoey();

    thread *pugzer1;
    thread *pugzer2;
    if (mOptions->usePugz) {
        pugzer1 = new thread(bind(&BarcodeToPositionMulti::pugzTask1, this));
        pugzer2 = new thread(bind(&BarcodeToPositionMulti::pugzTask2, this));
    }


    thread *getMbp;
    getMbp = new thread(bind(&BarcodeToPositionMulti::getMbpmap, this));


    getMbp->join();
#ifdef PRINT_INFO
    printf("processor %d get map thread done,cost %.4f\n", mOptions->myRank, GetTime() - t0);
#endif
    t0 = GetTime();

    thread *producer;
    producer = new thread(bind(&BarcodeToPositionMulti::producerTask, this));

    mOptions->dims1Size = mbpmap->GetDims1();


    Result **results = new Result *[mOptions->thread];
    BarcodeProcessor **barcodeProcessors = new BarcodeProcessor *[mOptions->thread];
    for (int t = 0; t < mOptions->thread; t++) {
        results[t] = new Result(mOptions, true);
//        results[t]->setBarcodeProcessor(mbpmap->GetHashNum(), mbpmap->GetHashHead(), mbpmap->GetHashMap(),
//                                        mbpmap->GetBloomFilter());
        results[t]->setBarcodeProcessorHashTableOneArrayWithBloomFilter(mbpmap->getHead(), mbpmap->getNext(),
                                                                        mbpmap->getPositionAll(),
                                                                        mbpmap->getBloomFilter());

    }
#ifdef PRINT_INFO
    printf("processor %d get results done,cost %.4f\n", mOptions->myRank, GetTime() - t0);
#endif
//    int threadsNowNumber = mOptions->thread - mOptions->pugzThread * 2;
    int threadsNowNumber = mOptions->thread;
    t0 = GetTime();
#ifdef PRINT_INFO
    printf("now start %d\n", threadsNowNumber);
#endif
    thread **threads = new thread *[mOptions->thread];
    for (int t = 0; t < threadsNowNumber; t++) {
        threads[t] = new thread(bind(&BarcodeToPositionMulti::consumerTask, this, results[t]));
    }


    thread *writerThread = NULL;
    thread *unMappedWriterThread = NULL;
    thread *mergeThread = NULL;
    if (mWriter) {
        writerThread = new thread(bind(&BarcodeToPositionMulti::writeTask, this, mWriter));
        if (mOptions->outGzSpilt == 0 && mOptions->numPro > 1 && mOptions->myRank == 0) {
            mergeThread = new thread(bind(&BarcodeToPositionMulti::mergeWrite, this));
        }
    }
    if (mUnmappedWriter) {
        unMappedWriterThread = new thread(bind(&BarcodeToPositionMulti::writeTask, this, mUnmappedWriter));
    }

    thread *pigzThread;
    if (mOptions->outGzSpilt) {
        pigzThread = new thread(bind(&BarcodeToPositionMulti::pigzWrite, this));
    } else {
        if (mOptions->usePigz && mOptions->myRank == 0) {
            pigzThread = new thread(bind(&BarcodeToPositionMulti::pigzWrite, this));
        }
    }
    if (mOptions->usePugz) {
        pugzer1->join();
        pugzer2->join();
    }
    double tt0 = GetTime();
#ifdef PRINT_INFO

    printf("now start another %d threads\n", mOptions->thread - threadsNowNumber);
#endif
//    for (int t = threadsNowNumber; t < mOptions->thread; t++) {
//        threads[t] = new thread(bind(&BarcodeToPositionMulti::consumerTask, this, results[t]));
//    }


    producer->join();

    for (int t = 0; t < mOptions->thread; t++) {
        threads[t]->join();
    }
#ifdef PRINT_INFO

    printf("another threads cost %.6f\n", GetTime() - tt0);
    printf("processor %d consumer cost %.4f\n", mOptions->myRank, GetTime() - t0);
#endif
    if (writerThread) {
        if (mOptions->outGzSpilt == 0 && mOptions->numPro > 1 && mOptions->myRank == 0) {
            mergeThread->join();
            mergeDone = 1;
        }
        writerThread->join();
        writerDone = 1;
#ifdef PRINT_INFO
        printf("processor %d writer done, cost %.4f\n", mOptions->myRank, GetTime() - t0);
#endif

        if (mOptions->outGzSpilt) {
            pigzThread->join();
#ifdef PRINT_INFO
            printf("processor %d pigz done, cost %.4f\n", mOptions->myRank, GetTime() - t0);
#endif
        } else {
            if (mOptions->usePigz && mOptions->myRank == 0) {
                pigzThread->join();
#ifdef PRINT_INFO
                printf("pigz done, cost %.4f\n", GetTime() - t0);
#endif
            }
        }
    }
    if (unMappedWriterThread)
        unMappedWriterThread->join();


    if (mOptions->verbose)
        loginfo("start to generate reports\n");

    //merge result
    vector<Result *> resultList;
    for (int t = 0; t < mOptions->thread; t++) {
        resultList.push_back(results[t]);
    }
    Result *finalResult = Result::merge(resultList);
#ifdef PRINT_INFO

    printf("###processor %d wait cost %f\n", mOptions->myRank, finalResult->GetCostWait());
    printf("###processor %d format cost %f\n", mOptions->myRank, finalResult->GetCostFormat());
    printf("###processor %d new cost %f\n", mOptions->myRank, finalResult->GetCostNew());
    printf("###processor %d pe cost %f\n", mOptions->myRank, finalResult->GetCostPe());
    printf("###processor %d all cost %f\n", mOptions->myRank, finalResult->GetCostAll());
#endif
    if (mOptions->numPro == 1) {
#ifdef PRINT_INFO

        printf("=======================print ans from process 0=========================\n");
#endif
        finalResult->print();
#ifdef PRINT_INFO

        printf("========================================================================\n");
#endif
    } else {
#ifdef PRINT_INFO

        printf("watind barrier\n");
#endif
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef PRINT_INFO

        printf("all merge done\n");
#endif
        if (mOptions->myRank == 0) {
#ifdef PRINT_INFO

            printf("=======================print ans from process 0=========================\n");
#endif
            vector<Result *> newResList;
            newResList.push_back(finalResult);
            for (int ii = 1; ii < mOptions->numPro; ii++) {
                Result *resultTmp = new Result(mOptions, true);
                resultTmp->setBarcodeProcessorHashTableOneArrayWithBloomFilter(mbpmap->getHead(), mbpmap->getNext(),
                                                                               mbpmap->getPositionAll(),
                                                                               mbpmap->getBloomFilter());
                MPI_Recv(&(resultTmp->mTotalRead), 1, MPI_LONG_LONG, ii, 1, mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mFxiedFilterRead), 1, MPI_LONG_LONG, ii, 1, mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mDupRead), 1, MPI_LONG_LONG, ii, 1, mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mLowQuaRead), 1, MPI_LONG_LONG, ii, 1, mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mWithoutPositionReads), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->overlapReadsWithMis), 1, MPI_LONG_LONG, ii, 1, mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->overlapReadsWithN), 1, MPI_LONG_LONG, ii, 1, mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->totalReads), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->mMapToSlideRead), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->overlapReads), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->overlapReadsWithMis), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->overlapReadsWithN), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->barcodeQ10), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->barcodeQ20), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->barcodeQ30), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->umiQ10), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->umiQ20), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->umiQ30), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->umiQ10FilterReads), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->umiNFilterReads), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->umiPloyAFilterReads), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->totQuery), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->filterQuery), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&(resultTmp->mBarcodeProcessor->queryYes), 1, MPI_LONG_LONG, ii, 1,
                         mOptions->communicator,
                         MPI_STATUS_IGNORE);
                newResList.push_back(resultTmp);
            }
            finalResult = Result::merge(newResList);
            finalResult->print();
#ifdef PRINT_INFO

            printf("========================================================================\n");
#endif
        } else {
            MPI_Send(&(finalResult->mTotalRead), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mFxiedFilterRead), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mDupRead), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mLowQuaRead), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mWithoutPositionReads), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->overlapReadsWithMis), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->overlapReadsWithN), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->totalReads), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->mMapToSlideRead), 1, MPI_LONG_LONG, 0, 1,
                     mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->overlapReads), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->overlapReadsWithMis), 1, MPI_LONG_LONG, 0, 1,
                     mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->overlapReadsWithN), 1, MPI_LONG_LONG, 0, 1,
                     mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->barcodeQ10), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->barcodeQ20), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->barcodeQ30), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->umiQ10), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->umiQ20), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->umiQ30), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->umiQ10FilterReads), 1, MPI_LONG_LONG, 0, 1,
                     mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->umiNFilterReads), 1, MPI_LONG_LONG, 0, 1,
                     mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->umiPloyAFilterReads), 1, MPI_LONG_LONG, 0, 1,
                     mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->totQuery), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->filterQuery), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
            MPI_Send(&(finalResult->mBarcodeProcessor->queryYes), 1, MPI_LONG_LONG, 0, 1, mOptions->communicator);
        }
    }


    cout << resetiosflags(ios::fixed) << setprecision(2);
    if (!mOptions->transBarcodeToPos.mappedDNBOutFile.empty()) {
        cout << "mapped_dnbs: " << finalResult->mBarcodeProcessor->mDNB.size() << endl;
        finalResult->dumpDNBs(mOptions->transBarcodeToPos.mappedDNBOutFile);
    }

    //clean up
    for (int t = 0; t < mOptions->thread; t++) {
        delete threads[t];
        threads[t] = NULL;
        delete results[t];
        results[t] = NULL;
    }

    delete[] threads;
    delete[] results;

    if (writerThread)
        delete writerThread;
    if (unMappedWriterThread)
        delete unMappedWriterThread;

    closeOutput();
#ifdef PRINT_INFO
    if (mOptions->myRank == 0)
        printf("final and delete cost %.4f\n", GetTime() - t0);
#endif
    return true;
}

void BarcodeToPositionMulti::initOutput() {
    mWriter = new WriterThread(mOptions->out, mOptions, mOptions->compression);
    if (!mOptions->transBarcodeToPos.unmappedOutFile.empty()) {
        mUnmappedWriter = new WriterThread(mOptions->transBarcodeToPos.unmappedOutFile, mOptions->compression);
    }
}

void BarcodeToPositionMulti::closeOutput() {
    if (mWriter) {
        delete mWriter;
        mWriter = NULL;
    }
    if (mUnmappedWriter) {
        delete mUnmappedWriter;
        mUnmappedWriter = NULL;
    }
}

bool BarcodeToPositionMulti::processPairEnd(ReadPairPack *pack, Result *result) {
    string outstr;
    string unmappedOut;
    bool hasPosition;
    bool fixedFiltered;
    for (int p = 0; p < pack->count; p++) {
        result->mTotalRead++;
        ReadPair *pair = pack->data[p];
        Read *or1 = pair->mLeft;
        Read *or2 = pair->mRight;
        if (filterFixedSequence) {
            fixedFiltered = fixedFilter->filter(or1, or2, result);
            if (fixedFiltered) {
                delete pair;
                continue;
            }
        }
        hasPosition = result->mBarcodeProcessor->process(or1, or2);
//        hasPosition = 1;
        if (hasPosition) {
            outstr += or2->toString();
        } else if (mUnmappedWriter) {
            unmappedOut += or2->toString();
        }
        delete pair;
    }
    mOutputMtx.lock();
    if (mUnmappedWriter && !unmappedOut.empty()) {
        //write reads that can't be mapped to the slide
        char *udata = new char[unmappedOut.size()];
        memcpy(udata, unmappedOut.c_str(), unmappedOut.size());
        mUnmappedWriter->input(udata, unmappedOut.size());
    }
    if (mWriter && !outstr.empty()) {
        char *data = new char[outstr.size()];
        memcpy(data, outstr.c_str(), outstr.size());
        //TODO add pigz queue here
        mWriter->input(data, outstr.size());
    }
    mOutputMtx.unlock();
//    delete pack->data;
    delete pack;
    return true;
}

void BarcodeToPositionMulti::initPackRepositoey() {
    mRepo.packBuffer = new ChunkPair *[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPairPack *) * PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
}

void BarcodeToPositionMulti::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void BarcodeToPositionMulti::producePack(ChunkPair *pack) {
    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void BarcodeToPositionMulti::consumePack(Result *result) {
    double tt = GetTime();

    double t = GetTime();
    dsrc::fq::FastqDataChunk *chunk;
    ChunkPair *chunkpair;
    ReadPairPack *data = new ReadPairPack;
    ReadPack *leftPack = new ReadPack;
    ReadPack *rightPack = new ReadPack;
    mInputMutx.lock();
    int cnt = 0;
    while (mRepo.writePos <= mRepo.readPos) {
        usleep(100);
        cnt++;
        if (mProduceFinished) {
            mInputMutx.unlock();
            return;
        }
    }
//    printf("P wait C %d times\n", cnt);
    chunkpair = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;
    mInputMutx.unlock();
    result->costWait += GetTime() - t;


    t = GetTime();
    leftPack->count = dsrc::fq::chunkFormat(chunkpair->leftpart, leftPack->data, true);
    rightPack->count = dsrc::fq::chunkFormat(chunkpair->rightpart, rightPack->data, true);
    result->costFormat += GetTime() - t;


    t = GetTime();
    data->count = leftPack->count < rightPack->count ? leftPack->count : rightPack->count;
    for (int i = 0; i < data->count; ++i) {
        data->data.push_back(new ReadPair(leftPack->data[i], rightPack->data[i]));
    }
    pairReader->fastqPool_left->Release(chunkpair->leftpart);
    pairReader->fastqPool_right->Release(chunkpair->rightpart);
    result->costNew += GetTime() - t;


    t = GetTime();
    processPairEnd(data, result);
    result->costPE += GetTime() - t;


    delete leftPack;
    delete rightPack;
    result->costAll += GetTime() - tt;
}


void BarcodeToPositionMulti::pugzTask1() {
#ifdef PRINT_INFO

    printf("now use pugz0 to decompress(%d threads)\n", mOptions->pugzThread);
#endif
    double t0 = GetTime();
    struct file_stream in;
    stat_t stbuf;
    int ret;
    const byte *in_p;

    ret = xopen_for_read(mOptions->transBarcodeToPos.in1.c_str(), true, &in);
    if (ret != 0) {
        printf("gg on xopen_for_read\n");
        exit(0);
    }

    ret = stat_file(&in, &stbuf, true);
    if (ret != 0) {
        printf("gg on stat_file\n");
        exit(0);
    }
    /* TODO: need a streaming-friendly solution */
    ret = map_file_contents(&in, size_t(stbuf.st_size));

    if (ret != 0) {
        printf("gg on map_file_contents\n");
        exit(0);
    }

    in_p = static_cast<const byte *>(in.mmap_mem);
    OutputConsumer output{};


    output.P = pugzQueue1;
    output.num = 1;
    output.pDone = &producerDone;
    ConsumerSync sync{};
    libdeflate_gzip_decompress(in_p, in.mmap_size, mOptions->pugzThread, output, &sync);

    pugz1Done = 1;
#ifdef PRINT_INFO

    std::cout << "pugz0 done, cost " << GetTime() - t0 << std::endl;
#endif
//    xclose(&in);
}

void BarcodeToPositionMulti::pugzTask2() {
#ifdef PRINT_INFO

    printf("now use pugz1 to decompress(%d threads)\n", mOptions->pugzThread);
#endif
    double t0 = GetTime();
    struct file_stream in;
    stat_t stbuf;
    int ret;
    const byte *in_p;

    ret = xopen_for_read(mOptions->transBarcodeToPos.in2.c_str(), true, &in);
    if (ret != 0) {
        printf("gg on xopen_for_read\n");
        exit(0);
    }

    ret = stat_file(&in, &stbuf, true);
    if (ret != 0) {
        printf("gg on stat_file\n");
        exit(0);
    }
    /* TODO: need a streaming-friendly solution */
    ret = map_file_contents(&in, size_t(stbuf.st_size));

    if (ret != 0) {
        printf("gg on map_file_contents\n");
        exit(0);
    }

    in_p = static_cast<const byte *>(in.mmap_mem);

    OutputConsumer output{};
    output.P = pugzQueue2;
    output.num = 2;
    output.pDone = &producerDone;
    ConsumerSync sync{};
    libdeflate_gzip_decompress(in_p, in.mmap_size, mOptions->pugzThread, output, &sync);
    pugz2Done = 1;
#ifdef PRINT_INFO

    std::cout << "pugz1 done, cost " << GetTime() - t0 << std::endl;
#endif
//    xclose(&in);
}

#define whoNumber 3

void BarcodeToPositionMulti::producerTask() {
    double t0 = GetTime();
    if (mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    pairReader = new FastqChunkReaderPair(mOptions->transBarcodeToPos.in1, mOptions->transBarcodeToPos.in2, true, 0, 0);

    ChunkPair *chunk_pair;
    pair<char *, int> last1;
    pair<char *, int> last2;
    int cnt = 0;
    int whoTurn = 0;
    int mps[whoNumber];
    if (mOptions->numPro == 1) {
        for (int i = 0; i < whoNumber; i++) {
            mps[i] = 0;
        }
    } else {
        mps[0] = 0;
        mps[1] = 1;
        mps[2] = 1;
//        mps[3] = 1 ;
//        mps[4] = 1 ;
//        mps[5] = 1 ;
//        mps[6] = 1 ;

    }
#ifdef PRINT_INFO

    printf("mps\n");
    for (int i = 0; i < whoNumber; i++)
        printf("%d ", mps[i]);
    printf("\n");
#endif
    long long p1Sum = 0;
    long long p2Sum = 0;
    if (mOptions->usePugz) {
        //TODO is this enough?
        last1.first = new char[1 << 20];
        last1.second = 0;
        last2.first = new char[1 << 20];
        last2.second = 0;


        while ((chunk_pair = pairReader->readNextChunkPair(pugzQueue1, pugzQueue2, &pugz1Done, &pugz2Done,
                                                           last1, last2)) != NULL) {
            //cerr << (char*)chunk_pair->leftpart->data.Pointer();
            if (mOptions->verbose)
                loginfo("producer read one chunk");
//            printf("read %d chunk done, size is %lld %lld\n", cnt++,
//                   chunk_pair->leftpart->size, chunk_pair->rightpart->size);
//            cout << "read " << (cnt++) << " chunk done, size is " << chunk_pair->leftpart->size << " "
//                 << chunk_pair->rightpart->size << endl;
            p1Sum += chunk_pair->leftpart->size;
            p2Sum += chunk_pair->rightpart->size;

            if (mps[whoTurn] == mOptions->myRank) {
                producePack(chunk_pair);
                cnt++;
            } else {
                pairReader->fastqPool_left->Release(chunk_pair->leftpart);
                pairReader->fastqPool_right->Release(chunk_pair->rightpart);
            }
            whoTurn++;
            whoTurn %= whoNumber;
            while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
//                printf("producer wait consumer\n");
//                cout << "producer wait consumer" << endl;
                slept++;
                usleep(100);
            }

        }

    } else {
        while ((chunk_pair = pairReader->readNextChunkPair()) != NULL) {
            //cerr << (char*)chunk_pair->leftpart->data.Pointer();
            if (mOptions->verbose)
                loginfo("producer read one chunk");
//            printf("read %d chunk done, size is %lld %lld\n", cnt++,
//                   chunk_pair->leftpart->size, chunk_pair->rightpart->size);
            p1Sum += chunk_pair->leftpart->size;
            p2Sum += chunk_pair->rightpart->size;
            if (mps[whoTurn] == mOptions->myRank) {
                producePack(chunk_pair);
                cnt++;
            } else {
                pairReader->fastqPool_left->Release(chunk_pair->leftpart);
                pairReader->fastqPool_right->Release(chunk_pair->rightpart);
            }
            whoTurn++;
            whoTurn %= whoNumber;
            while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
#ifdef PRINT_INFO
                printf("producer waiting...\n");
#endif
                slept++;
                usleep(100);
            }

        }
    }
#ifdef PRINT_INFO
    printf("processor %d  producer get1 %lld data\n", mOptions->myRank, p1Sum);
    printf("processor %d  producer get2 %lld data\n", mOptions->myRank, p2Sum);
#endif
    mProduceFinished = true;
    if (mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();
#ifdef PRINT_INFO
    printf("processor %d producer get %d chunk done, cost %.5f\n", mOptions->myRank, cnt, GetTime() - t0);
#endif
    producerDone = 1;
}


void BarcodeToPositionMulti::consumerTask(Result *result) {
    //cout << "in consumerTask " << endl;
    while (true) {
        while (mRepo.writePos <= mRepo.readPos) {
            if (mProduceFinished)
                break;
//            printf("consumer waiting...\n");
            usleep(100);
        }
        //unique_lock<mutex> lock(mRepo.readCounterMtx);
        if (mProduceFinished && mRepo.writePos == mRepo.readPos) {
            mFinishedThreads++;
            if (mOptions->verbose) {
                string msg = "finished " + to_string(mFinishedThreads) + " threads. Data processing completed.";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if (mProduceFinished) {
            if (mOptions->verbose) {
                string msg = "thread is processing the " + to_string(mRepo.readPos) + "/" + to_string(mRepo.writePos) +
                             " pack";
                loginfo(msg);
            }
            consumePack(result);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(result);
        }
    }


    if (mFinishedThreads == mOptions->thread) {
        if (mWriter)
            mWriter->setInputCompleted();
        if (mUnmappedWriter)
            mUnmappedWriter->setInputCompleted();
    }

    if (mOptions->verbose) {
        string msg = "finished one thread";
        loginfo(msg);
    }
}


void BarcodeToPositionMulti::writeTask(WriterThread *config) {
    if (mOptions->outGzSpilt) {
        while (true) {
            if (config->isCompleted()) {
                config->output(pigzQueue);
                break;
            }
            config->output(pigzQueue);
        }
#ifdef PRINT_INFO

        printf("processor %d wSum is %ld\n", mOptions->myRank, config->GetWSum());
#endif
    } else {
        if (mOptions->usePigz) {
            if (mOptions->myRank == 0) {
                while (true) {
                    if (config->isCompleted() && mergeDone) {
                        config->output(pigzQueue);
                        break;
                    }
                    config->output(pigzQueue);
                }
#ifdef PRINT_INFO

                printf("processor %d wSum is %lld\n", mOptions->myRank, config->GetWSum());
                printf("processor %d cSum is %d\n", mOptions->myRank, config->GetCSum());
#endif
            } else {
                while (true) {
                    if (config->isCompleted()) {
                        config->output(mOptions->communicator);
                        break;
                    }
                    config->output(mOptions->communicator);
                }
                int tag = -1;
#ifdef PRINT_INFO

                printf("processor 1 send data done, now send -1\n");
#endif
                MPI_Send(&(tag), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//                MPI_Barrier(MPI_COMM_WORLD);
#ifdef PRINT_INFO

                printf("processor 1 send -1 done\n");
#endif
            }
        } else {
            if (mOptions->myRank == 0) {
                while (true) {
                    if (config->isCompleted() && mergeDone) {
                        config->output();
                        break;
                    }
                    config->output();
                }
#ifdef PRINT_INFO

                printf("processor %d wSum is %lld\n", mOptions->myRank, config->GetWSum());
                printf("processor %d cSum is %d\n", mOptions->myRank, config->GetCSum());
#endif
            } else {
                while (true) {
                    if (config->isCompleted()) {
                        config->output(mOptions->communicator);
                        break;
                    }
                    config->output(mOptions->communicator);
                }
                int tag = -1;
#ifdef PRINT_INFO

                printf("processor 1 send data done, now send -1\n");
#endif
                MPI_Send(&(tag), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//                MPI_Barrier(MPI_COMM_WORLD);
#ifdef PRINT_INFO

                printf("processor 1 send -1 done\n");
#endif
            }
        }
    }

    if (mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
