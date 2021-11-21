/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
/*
  This file is modified by SDU HPC lab for RabbitQC project.

  Last modified: JULY2019 
*/
#ifndef H_FASTQSTREAM
#define H_FASTQSTREAM

#include <unistd.h>
#include <cstring>
#include "Globals.h"

#include "Fastq.h"
#include "Buffer.h"
#include "zlib/zlib.h"
#include "util.h"
#include "readerwriterqueue.h"
#include "atomicops.h"

#if defined (_WIN32)
#   define _CRT_SECURE_NO_WARNINGS
#   pragma warning(disable : 4996) // D_SCL_SECURE
#   pragma warning(disable : 4244) // conversion uint64 to uint32
//# pragma warning(disable : 4267)
#   define FOPEN    fopen
#   define FSEEK    _fseeki64
#   define FTELL    _ftelli64
#   define FCLOSE   fclose
#elif __APPLE__ // Apple by default suport 64 bit file operations (Darwin 10.5+)
#   define FOPEN    fopen
#   define FSEEK    fseek
#   define FTELL    ftell
#   define FCLOSE   fclose
#else
#   if !defined(_LARGEFILE_SOURCE)
#       define _LARGEFILE_SOURCE
#       if !defined(_LARGEFILE64_SOURCE)
#           define _LARGEFILE64_SOURCE
#       endif
#   endif
#   if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#       undef _FILE_OFFSET_BITS
#   endif
#   if !defined(_FILE_OFFSET_BITS)
#       define _FILE_OFFSET_BITS 64
#   endif
#   define FOPEN    fopen64
#   define FSEEK    fseeko64
#   define FTELL    ftello64
#   define FCLOSE   fclose
#endif
static const uint32 SwapBufferSize = 1 << 22;
static const uint32 tmpSwapBufferSize = 1 << 20;

namespace dsrc {

    namespace fq {

        class FastqFileReader {
        private:

        public:
            FastqFileReader(const std::string &fileName_)
                    : swapBuffer(SwapBufferSize), bufferSize(0), eof(false), usesCrlf(false), isZipped(false) {
                if (ends_with(fileName_, ".gz")) {
                    mZipFile = gzopen(fileName_.c_str(), "r");
                    isZipped = true;
                    gzrewind(mZipFile);

                } else {
                    mFile = FOPEN(fileName_.c_str(), "rb");
                    if (mFile == NULL) {
                        throw DsrcException(("Can not open file to read: " +
                                             fileName_).c_str()); //--------------need to change----------//
                    }
                }


            }

            ~FastqFileReader() {
                if (mFile != NULL || mZipFile != NULL)
                    Close();
                delete mFile;
                //delete mZipFile;
            }

            bool Eof() const {
                return eof;
            }

            bool ReadNextChunk(FastqDataChunk *chunk_);

            bool ReadNextPairedChunk(FastqDataChunk *chunk_);

            void Close() {
                if (mFile != NULL)
                    FCLOSE(mFile);
                if (mZipFile != NULL && isZipped) {
                    gzclose(mZipFile);
                    mZipFile = NULL;
                }

                mFile = NULL;
            }

            int64 Read(byte *memory_, uint64 size_) {
                if (isZipped) {
                    int64 n = gzread(mZipFile, memory_, size_);
                    if (n == -1)
                        cerr << "Error to read gzip file" << endl;
                    return n;
                } else {
                    int64 n = fread(memory_, 1, size_, mFile);
                    return n;
                }
            }

            int64 Read(byte *buf, uint64 len, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *Q,
                       atomic_int *done, pair<char *, int> &L, int num) {
                pair<char *, int> now;
                int64 ret;
                int64 got = 0;
                if (L.second > 0) {
                    if (L.second >= len) {
                        memcpy(buf, L.first, len);
                        char *tmp = new char[L.second - len];
                        memcpy(tmp, L.first + len, L.second - len);
                        memcpy(L.first, tmp, L.second - len);
                        delete[]tmp;
                        L.second = L.second - len;
                        ret = len;
                        buf += ret;
                        len -= ret;
                        got += ret;
                        return got;
                    } else {
                        memcpy(buf, L.first, L.second);
                        ret = L.second;
                        L.second = 0;
                        buf += ret;
                        len -= ret;
                        got += ret;
                    }
                }
                bool overWhile = false;
                while (len > 0) {
                    while (Q->try_dequeue(now) == 0) {
                        if (Q->size_approx() == 0 && *done == 1) {
                            ret = 0;
                            overWhile = true;
                            break;
                        }
                        usleep(100);
                    }
                    if (overWhile) {
                        ret = 0;
                        break;
                    }
//                    printf("get some data %d\n", now.second);
                    if (now.second <= len) {
                        memcpy(buf, now.first, now.second);
                        delete[] now.first;
                        ret = now.second;
                    } else {
                        int move_last = now.second - len;
                        memcpy(buf, now.first, len);
                        memcpy(L.first, now.first + len, move_last);
                        L.second = move_last;
                        delete[] now.first;
                        ret = len;
                    }

                    buf += ret;
                    len -= ret;
                    got += ret;
                }

                return got;
//
//                int64 resSum = 0;
//                int64 numFor = 0;
//                //TODO lastInfo.size > size_
//                if (lastInfo.second > size_) {
//                    printf("lastInfo.second>size_\n");
//                    exit(0);
//                }
//                if (lastInfo.second != 0) {
//                    memcpy(memory_, lastInfo.first, lastInfo.second);
//                    resSum += lastInfo.second;
//                    memory_ += lastInfo.second;
//                }
//                int cntt = 0;
//
//                while (resSum < size_) {
//
////                    if (numFor > 5000) {
////                        printf("num for > 5000 ,break\n");
////                        break;
////                    }
//
//
//                    std::pair<char *, int> now;
//                    if (q->try_dequeue(now)) {
//
//                        cntt = 0;
////                        cout << "get one small chunk, size is " << now.second << "queue size " << q->size_approx()
////                             << endl;
//                        if (resSum + now.second <= size_) {
////                            printf("read all small chunk to memory_...\n");
//                            memcpy(memory_, now.first, now.second);
//                            resSum += now.second;
//                            memory_ += now.second;
////                            printf("read after res sum is %lld\n\n", resSum);
//                            delete[] now.first;
//                        } else {
////                            printf("read part of small chunk to memory_...\n");
//                            int canCpoy = size_ - resSum;
//                            memcpy(memory_, now.first, canCpoy);
//                            resSum += canCpoy;
//                            memory_ += canCpoy;
//                            lastInfo.second = now.second - canCpoy;
//                            memcpy(lastInfo.first, now.first + canCpoy, lastInfo.second);
////                            printf("read after res sum is %lld\n\n", resSum);
//
//                            delete[] now.first;
//                            break;
//                        }
//                    } else if (q->size_approx() == 0) {
//                        if (done == 1)break;
//                        numFor++;
////                        printf("prudocer wait pugz %d\n", cntt++);
////                        cout << "prudocer " << num << " wait pugz " << cntt++ << " q->size_approx() "
////                             << q->size_approx() << endl;
//                        usleep(100);
//                    }
//                }
////                printf("ok of one read, now res sum is %lld, need to read is %lld, leave %d to read next time\n",
////                       resSum, size_, lastInfo.second);
////                printf("===============================\n");
//                return resSum;
            }

        private:
            core::Buffer swapBuffer;
            uint64 bufferSize;
            bool eof;
            bool usesCrlf;
            bool isZipped;
            FILE *mFile;
            gzFile mZipFile;

            uint64 lastOneReadPos;
            uint64 lastTwoReadPos;

            uint64 GetNextRecordPos(uchar *data_, uint64 pos_, const uint64 size_);

            uint64 GetPreviousRecordPos(uchar *data_, uint64 pos_, const uint64 size_);

            void SkipToEol(uchar *data_, uint64 &pos_, const uint64 size_) {
                //cout << "SkipToEol " << pos_ << " " << size_ << endl;
                ASSERT(pos_ < size_);

                while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
                    ++pos_;

                if (data_[pos_] == '\r' && pos_ < size_) {
                    if (data_[pos_ + 1] == '\n') {
                        usesCrlf = true;
                        ++pos_;
                    }
                }
            }

            //跳转到行首
            void SkipToSol(uchar *data_, uint64 &pos_, const uint64 size_) {
                //std::cout<<"SkipToSol:"<<data_[pos_]<<std::endl;
                ASSERT(pos_ < size_);
                if (data_[pos_] == '\n') {
                    --pos_;
                }
                if (data_[pos_] == '\r') {
                    usesCrlf = true;
                    pos_--;
                }
                //找到换行符
                while (data_[pos_] != '\n' && data_[pos_] != '\r') {
                    //std::cout<<"pos_;"<<pos_<<std::endl;
                    --pos_;
                }
                if (data_[pos_] == '\n') {
                    --pos_;
                }
                if (data_[pos_] == '\r') {
                    usesCrlf = true;
                    pos_--;
                }

            }


        };


    } // namespace fq

} // namespace dsrc

#endif // H_FASTQSTREAM
