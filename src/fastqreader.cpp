#include "fastqreader.h"
#include "util.h"
#include <string.h>
#include "FastqStream.h"


#define FQ_BUF_SIZE (1<<20)

FastqReader::FastqReader(string filename, bool hasQuality, bool phred64) {
    mFilename = filename;
    mZipFile = NULL;
    mZipped = false;
    mFile = NULL;
    mStdinMode = false;
    mPhred64 = phred64;
    mHasQuality = hasQuality;
    mBuf = new char[FQ_BUF_SIZE];
    mBufDataLen = 0;
    mBufUsedLen = 0;
    mHasNoLineBreakAtEnd = false;
    init();
}

FastqReader::~FastqReader() {
    close();
    delete mBuf;
}

bool FastqReader::hasNoLineBreakAtEnd() {
    return mHasNoLineBreakAtEnd;
}

void FastqReader::readToBuf() {
    if (mZipped) {
        mBufDataLen = gzread(mZipFile, mBuf, FQ_BUF_SIZE);
        if (mBufDataLen == -1) {
            error_exit("Error to read gzip file: " + mFilename);
            //cerr << "Error to read gzip file" << endl;
        }
    } else {
        mBufDataLen = fread(mBuf, 1, FQ_BUF_SIZE, mFile);
    }
    mBufUsedLen = 0;

    if (mBufDataLen < FQ_BUF_SIZE) {
        if (mBuf[mBufDataLen - 1] != '\n')
            mHasNoLineBreakAtEnd = true;
    }
}

void FastqReader::init() {
    if (ends_with(mFilename, ".gz")) {
        mZipFile = gzopen(mFilename.c_str(), "r");
        mZipped = true;
        gzrewind(mZipFile);
    } else {
        if (mFilename == "/dev/stdin") {
            mFile = stdin;
        } else
            mFile = fopen(mFilename.c_str(), "rb");
        if (mFile == NULL) {
            error_exit("Failed to open file: " + mFilename);
        }
        mZipped = false;
    }
    readToBuf();
}

void FastqReader::getBytes(size_t &bytesRead, size_t &bytesTotal) {
    if (mZipped) {
        bytesRead = gzoffset(mZipFile);
    } else {
        bytesRead = ftell(mFile);//mFile.tellg();
    }

    // use another ifstream to not affect current reader
    ifstream is(mFilename);
    is.seekg(0, is.end);
    bytesTotal = is.tellg();
}

void FastqReader::clearLineBreaks(char *line) {

    // trim \n, \r or \r\n in the tail
    int readed = strlen(line);
    if (readed >= 2) {
        if (line[readed - 1] == '\n' || line[readed - 1] == '\r') {
            line[readed - 1] = '\0';
            if (line[readed - 2] == '\r')
                line[readed - 2] = '\0';
        }
    }
}

string FastqReader::getLine() {
    static int c = 0;
    c++;
    int copied = 0;

    int start = mBufUsedLen;
    int end = start;

    while (end < mBufDataLen) {
        if (mBuf[end] != '\r' && mBuf[end] != '\n')
            end++;
        else
            break;
    }

    // this line well contained in this buf, or this is the last buf
    if (end < mBufDataLen || mBufDataLen < FQ_BUF_SIZE) {
        int len = end - start;
        string line(mBuf + start, len);

        // skip \n or \r
        end++;
        // handle \r\n
        if (end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n')
            end++;

        mBufUsedLen = end;

        return line;
    }

    // this line is not contained in this buf, we need to read new buf
    string str(mBuf + start, mBufDataLen - start);

    while (true) {
        readToBuf();
        start = 0;
        end = 0;
        while (end < mBufDataLen) {
            if (mBuf[end] != '\r' && mBuf[end] != '\n')
                end++;
            else
                break;
        }
        // this line well contained in this buf, we need to read new buf
        if (end < mBufDataLen || mBufDataLen < FQ_BUF_SIZE) {
            int len = end - start;
            str.append(mBuf + start, len);

            // skip \n or \r
            end++;
            // handle \r\n
            if (end < mBufDataLen - 1 && mBuf[end] == '\n')
                end++;

            mBufUsedLen = end;
            return str;
        }
        // even this new buf is not enough, although impossible
        str.append(mBuf + start, mBufDataLen);
    }

    return string();
}

bool FastqReader::eof() {
    if (mZipped) {
        return gzeof(mZipFile);
    } else {
        return feof(mFile);//mFile.eof();
    }
}

Read *FastqReader::read() {
    if (mZipped) {
        if (mZipFile == NULL)
            return NULL;
    }

    if (mBufUsedLen >= mBufDataLen && eof()) {
        return NULL;
    }

    string name = getLine();
    // name should start with @
    while ((name.empty() && !(mBufUsedLen >= mBufDataLen && eof())) || (!name.empty() && name[0] != '@')) {
        name = getLine();
    }

    if (name.empty())
        return NULL;

    string sequence = getLine();
    string strand = getLine();

    // WAR for FQ with no quality
    if (!mHasQuality) {
        string quality = string(sequence.length(), 'K');
        return new Read(name, sequence, strand, quality, mPhred64);
    } else {
        string quality = getLine();
        if (quality.length() != sequence.length()) {
            cerr << "ERROR: sequence and quality have different length:" << endl;
            cerr << name << endl;
            cerr << sequence << endl;
            cerr << strand << endl;
            cerr << quality << endl;
            return NULL;
        }
        return new Read(name, sequence, strand, quality, mPhred64);
    }

    return NULL;
}

void FastqReader::close() {
    if (mZipped) {
        if (mZipFile) {
            gzclose(mZipFile);
            mZipFile = NULL;
        }
    } else {
        if (mFile) {
            fclose(mFile);//mFile.close();
            mFile = NULL;
        }
    }
}

bool FastqReader::isZipFastq(string filename) {
    if (ends_with(filename, ".fastq.gz"))
        return true;
    else if (ends_with(filename, ".fq.gz"))
        return true;
    else if (ends_with(filename, ".fasta.gz"))
        return true;
    else if (ends_with(filename, ".fa.gz"))
        return true;
    else
        return false;
}

bool FastqReader::isFastq(string filename) {
    if (ends_with(filename, ".fastq"))
        return true;
    else if (ends_with(filename, ".fq"))
        return true;
    else if (ends_with(filename, ".fasta"))
        return true;
    else if (ends_with(filename, ".fa"))
        return true;
    else
        return false;
}

bool FastqReader::isZipped() {
    return mZipped;
}

bool FastqReader::test() {
    FastqReader reader1("testdata/R1.fq");
    FastqReader reader2("testdata/R1.fq.gz");
    Read *r1 = NULL;
    Read *r2 = NULL;
    while (true) {
        r1 = reader1.read();
        r2 = reader2.read();
        if (r1 == NULL || r2 == NULL)
            break;
        if (r1->mSeq.mStr != r2->mSeq.mStr) {
            return false;
        }
        delete r1;
        delete r2;
    }
    return true;
}

FastqReaderPair::FastqReaderPair(FastqReader *left, FastqReader *right) {
    mLeft = left;
    mRight = right;
}

FastqReaderPair::FastqReaderPair(string leftName, string rightName, bool hasQuality, bool phred64, bool interleaved) {
    mInterleaved = interleaved;
    mLeft = new FastqReader(leftName, hasQuality, phred64);
    if (mInterleaved)
        mRight = NULL;
    else
        mRight = new FastqReader(rightName, hasQuality, phred64);
}

FastqReaderPair::~FastqReaderPair() {
    if (mLeft) {
        delete mLeft;
        mLeft = NULL;
    }
    if (mRight) {
        delete mRight;
        mRight = NULL;
    }
}

ReadPair *FastqReaderPair::read() {
    Read *l = mLeft->read();
    Read *r = NULL;
    if (mInterleaved)
        r = mLeft->read();
    else
        r = mRight->read();
    if (!l || !r) {
        return NULL;
    } else {
        return new ReadPair(l, r);
    }
}

//----------------------
FastqChunkReaderPair::FastqChunkReaderPair(dsrc::fq::FastqReader *left, dsrc::fq::FastqReader *right)
    : swapBuffer_left(SwapBufferSize), swapBuffer_right(SwapBufferSize), bufferSize_left(0), bufferSize_right(0), eof(false),
    usesCrlf(false) {
        mLeft = left;
        mRight = right;
    }

FastqChunkReaderPair::FastqChunkReaderPair(string leftName, string rightName, bool hasQuality, bool phred64,
        bool interleaved)
    : swapBuffer_left(SwapBufferSize), swapBuffer_right(SwapBufferSize), bufferSize_left(0), bufferSize_right(0), eof(false),
    usesCrlf(false) {
        mInterleaved = interleaved;
        fastqPool_left = new dsrc::fq::FastqDataPool(128, SwapBufferSize);
        fileReader_left = new dsrc::fq::FastqFileReader(leftName);
        mLeft = new dsrc::fq::FastqReader(*fileReader_left, *fastqPool_left);
        if (mInterleaved) {
            mRight = NULL;
            fastqPool_right = NULL;
            fileReader_right = NULL;
        } else {
            fastqPool_right = new dsrc::fq::FastqDataPool(128, SwapBufferSize);
            fileReader_right = new dsrc::fq::FastqFileReader(rightName);
            mRight = new dsrc::fq::FastqReader(*fileReader_right, *fastqPool_left);
        }
    }

FastqChunkReaderPair::~FastqChunkReaderPair() {
    if (mLeft) {
        delete fastqPool_left;
        delete fileReader_left;
        delete mLeft;
        fastqPool_left = NULL;
        fileReader_left = NULL;
        mLeft = NULL;
    }
    if (mRight) {
        delete fastqPool_right;
        delete fileReader_right;
        delete mRight;
        mRight = NULL;
    }
}

ChunkPair *FastqChunkReaderPair::readNextChunkPair() {
    if (mInterleaved) {
        //interleaved code here
        return NULL;

    } else {
        //not interleaved
        return readNextChunkPair_interleaved();
    }
}

ChunkPair *FastqChunkReaderPair::readNextChunkPair(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q1,
        moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q2,
        atomic_int *d1, atomic_int *d2,
        pair<char *, int> &last1, pair<char *, int> &last2) {
    if (mInterleaved) {
        //interleaved code here
        return NULL;

    } else {
        //not interleaved
        return readNextChunkPair_interleaved(q1, q2, d1, d2, last1, last2);
    }
}
//ChunkPair* FastqChunkReaderPair::readNextChunkPair(){
//ChunkPair* pair = new ChunkPair;
//dsrc::fq::FastqDataChunk* leftPart = NULL;
//fastqPool_left->Acquire(leftPart);
//
//dsrc::fq::FastqDataChunk* rightPart = NULL;
//fastqPool_right->Acquire(rightPart);
//
//if(fileReader_left -> ReadNextChunk(leftPart)){
////cerr << "left reader read one chunk size: " << leftPart->size << endl;
//pair->leftpart = leftPart;
////cerr << (char*) pair->leftpart->data.Pointer();
//}else{
//fastqPool_left->Release(leftPart);
//return NULL;
//}
//if(mInterleaved){
//if(fileReader_left -> ReadNextChunk(rightPart)){
//pair->rightpart = rightPart;
//}else{
//fastqPool_right->Release(rightPart);
//pair->rightpart = NULL;
//}
//}else{
//if(fileReader_right -> ReadNextChunk(rightPart)){
////cerr << "right reader read one chunk size: " << rightPart->size << endl;
//pair->rightpart = rightPart;
//}else{
//fastqPool_right->Release(rightPart);
//pair->rightpart = NULL;
//}
//}
//
//return pair;
//
//}
//uint64 count_line(dsrc::uchar* contenx, uint64 read_bytes){
//	__m128i current_v, mask_v, tmp_v;
//	uint8_t slitn = '\n';
//	__m128i slitn_v = _mm_set1_epi8(slitn);
//	__m128i one_v = _mm_set1_epi8(1);
//	__m128i zero_v = _mm_setzero_si128();
//	__m128i result_v = _mm_setzero_si128();
//	uint64 result = 0;
//	for(uint64 i = 0; i < read_bytes; i+=16){
//		current_v = _mm_load_si128((__m128i*)(contenx+i));
//		if(i%(16*256) == 0){
//			for(int j = 0; j < 16; j++){
//				result += ((uint8_t*)&result_v)[j];
//			}
//			result_v = _mm_setzero_si128();
//		}
//		mask_v = _mm_cmpeq_epi8(current_v, slitn_v);
//		tmp_v = _mm_blendv_epi8(zero_v, one_v, mask_v);
//		result_v = _mm_add_epi8(result_v, tmp_v);
//	}
//	for(int i = 0; i < 16; i++){
//		//printf("result_v[%d] is: %d\n", i, ((uint8_t*)&result_v)[i]);
//		result += ((uint8_t*)&result_v)[i];
//	}
//	return result;
//}

int64 count_line(dsrc::uchar *contenx, int64 read_bytes) {
    int64 count_n = 0;
    for (int i = 0; i < read_bytes; ++i) {
        //printf("%c",contenx[i]);
        if (contenx[i] == '\n') count_n++;
    }
    return count_n;
}

void FastqChunkReaderPair::SkipToEol(dsrc::uchar *data_, uint64 &pos_, const uint64 size_) {
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

uint64 FastqChunkReaderPair::GetNextRecordPos(dsrc::uchar *data_, uint64 pos_, const uint64 size_) {
    SkipToEol(data_, pos_, size_);
    ++pos_;

    // find beginning of the next record
    while (data_[pos_] != '@') {
        SkipToEol(data_, pos_, size_);
        ++pos_;
    }
    uint64 pos0 = pos_;

    SkipToEol(data_, pos_, size_);
    ++pos_;

    if (data_[pos_] == '@')            // previous one was a quality field
        return pos_;

    SkipToEol(data_, pos_, size_);
    ++pos_;
    if (data_[pos_] != '+')
        std::cout << "core dump is " << data_[pos_] << std::endl;
    ASSERT(data_[pos_] == '+');    // pos0 was the start of tag
    return pos0;
}

//ChunkPair* FastqChunkReaderPair::readNextChunkPair_allinone(){
ChunkPair *FastqChunkReaderPair::readNextChunkPair_interleaved() {
    int eof1=false;
    int eof2=false;
    ChunkPair *pair = new ChunkPair;
    dsrc::fq::FastqDataChunk *leftPart = NULL;
    fastqPool_left->Acquire(leftPart);

    dsrc::fq::FastqDataChunk *rightPart = NULL;
    fastqPool_right->Acquire(rightPart);

    int64 left_line_count = 0;
    int64 right_line_count = 0;
    int64 chunkEnd_right = 0;
    int64 chunkEnd = 0;
    //------read left chunk------//
    if (eof) {
        leftPart->size = 0;
        rightPart->size = 0;
        //return false;
        fastqPool_left->Release(leftPart);
        fastqPool_right->Release(rightPart);
        return NULL;
    }

    //flush the data from previous incomplete chunk
    dsrc::uchar *data = leftPart->data.Pointer();
    uint64 cbufSize = leftPart->data.Size();
    leftPart->size = 0;
    int64 toRead;
    toRead = cbufSize - bufferSize_left;
    if (bufferSize_left > 0) {
        /*
        int cntt = 0;
        int pos = 0;
        printf("===========OOOOOOL========\n");
        while (pos < bufferSize_left && cntt < 4) {
            printf("%c", swapBuffer_left.Pointer()[pos]);
            if (swapBuffer_left.Pointer()[pos] == '\n')cntt++;
            pos++;

        }
        printf("==========================\n");
        */    
        std::copy(swapBuffer_left.Pointer(), swapBuffer_left.Pointer() + bufferSize_left, data);
        leftPart->size = bufferSize_left;
        bufferSize_left = 0;
    }
    int64 r;
    r = mLeft->Read(data + leftPart->size, toRead);
    if (r > 0) {
        if (r == toRead) {
            chunkEnd = cbufSize - tmpSwapBufferSize;
            chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
        } else {
            printf("read1 r!=toRead\n");
            //chunkEnd = r;
            leftPart->size += r - 1;
            if (usesCrlf)
                leftPart->size -= 1;
            eof1 = true;
            chunkEnd = leftPart->size+1;
            cbufSize = leftPart->size+1;

        }
    } else {
        printf("read1 r <= 0\n");
        eof1 = true;
        chunkEnd = leftPart->size+1;
        cbufSize = leftPart->size+1;
        if(leftPart->size==0){
            return NULL;
        }
    }
    //------read left chunk end------//

    //-----------------read right chunk---------------------//
    dsrc::uchar *data_right = rightPart->data.Pointer();
    uint64 cbufSize_right = rightPart->data.Size();
    rightPart->size = 0;
    toRead = cbufSize_right - bufferSize_right;
    if (bufferSize_right > 0) {
        /*
        int cntt = 0;
        int pos = 0;
        printf("===========OOOOOOR========\n");
        while (pos < bufferSize_right && cntt < 4) {
            printf("%c", swapBuffer_right.Pointer()[pos]);
            if (swapBuffer_right.Pointer()[pos] == '\n')cntt++;
            pos++;

        }
        printf("==========================\n");
        */ 
        std::copy(swapBuffer_right.Pointer(), swapBuffer_right.Pointer() + bufferSize_right, data_right);
        rightPart->size = bufferSize_right;
        bufferSize_right = 0;
    }
    r = mRight->Read(data_right + rightPart->size, toRead);
    if (r > 0) {
        if (r == toRead) {
            chunkEnd_right = cbufSize_right - tmpSwapBufferSize;
            chunkEnd_right = GetNextRecordPos(data_right, chunkEnd_right, cbufSize_right);
        } else {
            printf("read2 r!=toRead\n");
            //chunkEnd_right += r;
            rightPart->size += r - 1;
            if (usesCrlf)
                rightPart->size -= 1;
            eof2 = true;
            chunkEnd_right = rightPart->size+1;
            cbufSize_right = rightPart->size+1;

        }
    } else {
        printf("read2 r>=0\n");
        eof2 = true;
        chunkEnd_right = rightPart->size+1;
        cbufSize_right = rightPart->size+1;
        if(rightPart->size==0)
            return NULL;
    }
    if(eof1 && eof2)eof = true;
    if(eof1 || eof2)printf("eofs : %d %d %d\n",eof1,eof2,eof);

    //--------------read right chunk end---------------------//
    if (!eof) {
        left_line_count = count_line(data, chunkEnd);
        right_line_count = count_line(data_right, chunkEnd_right);
        int64 difference = left_line_count - right_line_count;

        if (difference > 0) {
            while (chunkEnd >=0) {
                if (data[chunkEnd] == '\n') {
                    difference--;
                    if (difference == -1) {
                        chunkEnd++;
                        break;
                    }
                }
                chunkEnd--;
            }

        } else if (difference < 0) {

            while (chunkEnd_right >= 0) {
                if (data_right[chunkEnd_right] == '\n') {
                    difference++;
                    if (difference == 1) {
                        chunkEnd_right++;
                        break;
                    }
                }
                chunkEnd_right--;
            }
        }

        leftPart->size = chunkEnd - 1;
        if (usesCrlf)
            leftPart->size -= 1;

        std::copy(data + chunkEnd, data + cbufSize, swapBuffer_left.Pointer());
        bufferSize_left = cbufSize - chunkEnd;

        rightPart->size = chunkEnd_right - 1;
        if (usesCrlf)
            rightPart->size -= 1;
        std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer_right.Pointer());
        bufferSize_right = cbufSize_right - chunkEnd_right;
        left_line_count = count_line(data, chunkEnd);
        right_line_count = count_line(data_right, chunkEnd_right);
        //int64 difference = right_line_count - left_line_count;
        difference = left_line_count - right_line_count;
        if(difference!=0)printf("GG, diff is %d\n",difference); 
    }
    pair->leftpart = leftPart;
    pair->rightpart = rightPart;
    return pair;

}

ChunkPair *
FastqChunkReaderPair::readNextChunkPair_interleaved(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q1,
        moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q2,
        atomic_int *d1, atomic_int *d2,
        pair<char *, int> &last1, pair<char *, int> &last2) {
    ChunkPair *pair = new ChunkPair;
    dsrc::fq::FastqDataChunk *leftPart = NULL;
    fastqPool_left->Acquire(leftPart);
    int eof1=false;
    int eof2=false;

    dsrc::fq::FastqDataChunk *rightPart = NULL;
    fastqPool_right->Acquire(rightPart);

    int64 left_line_count = 0;
    int64 right_line_count = 0;
    int64 chunkEnd_right = 0;
    int64 chunkEnd = 0;
    //------read left chunk------//
    if (eof) {
        leftPart->size = 0;
        rightPart->size = 0;
        //return false;
        fastqPool_left->Release(leftPart);
        fastqPool_right->Release(rightPart);
        return NULL;
    }

    //flush the data from previous incomplete chunk
    dsrc::uchar *data = leftPart->data.Pointer();
    uint64 cbufSize = leftPart->data.Size();
    leftPart->size = 0;
    int64 toRead;
    toRead = cbufSize - bufferSize_left;
    if (bufferSize_left > 0) {
        std::copy(swapBuffer_left.Pointer(), swapBuffer_left.Pointer() + bufferSize_left, data);
        leftPart->size = bufferSize_left;
        bufferSize_left = 0;
    }
    int64 r;
    r = mLeft->Read(data + leftPart->size, toRead, q1, d1, last1, 1);
    //    printf("now read once done, read %lld / %lld\n", r, toRead);
    if (r > 0) {
        if (r == toRead) {
            chunkEnd = cbufSize - tmpSwapBufferSize;
            chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
            //leftPart->size = chunkEnd - 1;
            //if(usesCrlf)
            //	leftPart->size -= 1;

            //std::copy(data + chunkEnd, data + cbufSize, swapBuffer_left.Pointer());
            //bufferSize_left = cbufSize - chunkEnd;
        } else {
            //chunkEnd = r;
            leftPart->size += r - 1;
            if (usesCrlf)
                leftPart->size -= 1;
            eof1 = true;
            chunkEnd = leftPart->size+1;
            cbufSize = leftPart->size+1;
        }
    } else {
        eof1 = true;
        chunkEnd = leftPart->size+1;
        cbufSize = leftPart->size+1;
        if(leftPart->size==0){
            return NULL;
        }
    }
    //------read left chunk end------//

    //-----------------read right chunk---------------------//
    dsrc::uchar *data_right = rightPart->data.Pointer();
    uint64 cbufSize_right = rightPart->data.Size();
    rightPart->size = 0;
    toRead = cbufSize_right - bufferSize_right;
    if (bufferSize_right > 0) {
        std::copy(swapBuffer_right.Pointer(), swapBuffer_right.Pointer() + bufferSize_right, data_right);
        rightPart->size = bufferSize_right;
        bufferSize_right = 0;
    }
    r = mRight->Read(data_right + rightPart->size, toRead, q2, d2, last2, 2);
    //    printf("now read once done, read %lld / %lld\n", r, toRead);

    if (r > 0) {
        if ( r == toRead) {
            chunkEnd_right = cbufSize_right - tmpSwapBufferSize;
            chunkEnd_right = GetNextRecordPos(data_right, chunkEnd_right, cbufSize_right);
            //leftPart->size = chunkEnd - 1;
            //if(usesCrlf)
            //	leftPart->size -= 1;

            //std::copy(data + chunkEnd, data + cbufSize, swapBuffer_left.Pointer());
            //bufferSize_left = cbufSize - chunkEnd;
        } else {
            //chunkEnd_right += r;
            rightPart->size += r - 1;
            if (usesCrlf)
                rightPart->size -= 1;
            eof2 = true;
            chunkEnd_right = rightPart->size+1;
            cbufSize_right = rightPart->size+1;

        }
    } else {
        eof2 = true;
        chunkEnd_right = rightPart->size+1;
        cbufSize_right = rightPart->size+1;
        if(rightPart->size==0){
            return NULL;
        }
    }
    //--------------read right chunk end---------------------//
    if(eof1 && eof2)eof = true;
    if(eof1 || eof2)printf("eofs : %d %d %d\n",eof1,eof2,eof);
    if (!eof) {
        left_line_count = count_line(data, chunkEnd);
        right_line_count = count_line(data_right, chunkEnd_right);
        int64 difference = left_line_count - right_line_count;

        if (difference > 0) {
            while (chunkEnd >= 0) {
                if (data[chunkEnd] == '\n') {
                    difference--;
                    if (difference == -1) {
                        chunkEnd++;
                        break;
                    }
                }
                chunkEnd--;
            }

        } else if (difference < 0) {

            while (chunkEnd_right >= 0) {
                if (data_right[chunkEnd_right] == '\n') {
                    difference++;
                    if (difference == 1) {
                        chunkEnd_right++;
                        break;
                    }
                }
                chunkEnd_right--;
            }
        }
        leftPart->size = chunkEnd - 1;
        if (usesCrlf)
            leftPart->size -= 1;

        std::copy(data + chunkEnd, data + cbufSize, swapBuffer_left.Pointer());
        bufferSize_left = cbufSize - chunkEnd;

        rightPart->size = chunkEnd_right - 1;
        if (usesCrlf)
            rightPart->size -= 1;
        std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer_right.Pointer());
        bufferSize_right = cbufSize_right - chunkEnd_right;
        left_line_count = count_line(data, chunkEnd);
        right_line_count = count_line(data_right, chunkEnd_right);
        //int64 difference = right_line_count - left_line_count;
        difference = left_line_count - right_line_count;
        if(difference!=0)printf("GG, diff is %d\n",difference); 
    }
    pair->leftpart = leftPart;
    pair->rightpart = rightPart;
    return pair;

}


