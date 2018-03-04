#include <stdio.h>
#include <memory.h>
#include <inttypes.h>
#include <vector>
 
// typedefs
typedef uint16_t    uint16;
typedef uint32_t    uint32;
typedef int32_t     int32;
 
//this struct is the minimal required header data for a wav file
struct SMinimalWaveFileHeader
{
    //the main chunk
    unsigned char m_chunkID[4];
    uint32        m_chunkSize;
    unsigned char m_format[4];
 
    //sub chunk 1 "fmt "
    unsigned char m_subChunk1ID[4];
    uint32        m_subChunk1Size;
    uint16        m_audioFormat;
    uint16        m_numChannels;
    uint32        m_sampleRate;
    uint32        m_byteRate;
    uint16        m_blockAlign;
    uint16        m_bitsPerSample;
 
    //sub chunk 2 "data"
    unsigned char m_subChunk2ID[4];
    uint32        m_subChunk2Size;
 
    //then comes the data!
};


inline void FloatToPCM(unsigned char *PCM, const float& in, size_t numBytes)
{
    uint32 data;
    if (in < 0.0f)
        data = uint32(in * float(0x80000000));
    else
        data = uint32(in * float(0x7fffffff));

    switch (numBytes)
    {
        case 4: PCM[3] = ((data >> 24) & 0xFF); PCM[2] = ((data >> 16) & 0xFF); PCM[1] = ((data >> 8) & 0xFF); PCM[0] = (data & 0xFF); break;
        case 3: PCM[2] = ((data >> 24) & 0xFF); PCM[1] = ((data >> 16) & 0xFF); PCM[0] = ((data >> 8) & 0xFF); break;
        case 2: PCM[1] = ((data >> 24) & 0xFF); PCM[0] = ((data >> 16) & 0xFF); break; 
        case 1: PCM[0] = ((data >> 24) & 0xFF); break;
    }
}

inline void PCMToFloat(float& out, const unsigned char *PCM, size_t numBytes)
{
    uint32 data = 0;
    switch (numBytes)
    {
        case 4: data = (uint32(PCM[3]) << 24) | (uint32(PCM[2]) << 16) | (uint32(PCM[1]) << 8) | uint32(PCM[0]); break;
        case 3: data = (uint32(PCM[2]) << 24) | (uint32(PCM[1]) << 16) | (uint32(PCM[0]) << 8); break;
        case 2: data = (uint32(PCM[1]) << 24) | (uint32(PCM[0]) << 16); break;
        case 1: data = (uint32(PCM[0]) << 24); break;
    }

    if (data & 0x80000000)
        out = float(int32(data)) / float(0x80000000);
    else
        out = float(data) / float(0x7fffffff);
}

// TODO: make this function take bytes per sample, 1, 2, 3, or 4!
bool WriteWaveFile(const char *fileName, std::vector<float>& dataFloat, uint16 numChannels, uint32 sampleRate)
{
    std::vector<int32> data;
    data.resize(dataFloat.size());
    for (size_t i = 0; i < dataFloat.size(); ++i)
        FloatToPCM((unsigned char*)&data[i], dataFloat[i], 4);

    uint32 dataSize = (uint32)(data.size() * sizeof(int32));
    uint16 bitsPerSample = sizeof(int32) * 8;
 
    //open the file if we can
    FILE *File = nullptr;
    fopen_s(&File, fileName, "w+b");
    if (!File)
        return false;
 
    SMinimalWaveFileHeader waveHeader;
 
    //fill out the main chunk
    memcpy(waveHeader.m_chunkID, "RIFF", 4);
    waveHeader.m_chunkSize = dataSize + 36;
    memcpy(waveHeader.m_format, "WAVE", 4);
 
    //fill out sub chunk 1 "fmt "
    memcpy(waveHeader.m_subChunk1ID, "fmt ", 4);
    waveHeader.m_subChunk1Size = 16;
    waveHeader.m_audioFormat = 1;
    waveHeader.m_numChannels = numChannels;
    waveHeader.m_sampleRate = sampleRate;
    waveHeader.m_byteRate = sampleRate * numChannels * bitsPerSample / 8;
    waveHeader.m_blockAlign = numChannels * bitsPerSample / 8;
    waveHeader.m_bitsPerSample = bitsPerSample;
 
    //fill out sub chunk 2 "data"
    memcpy(waveHeader.m_subChunk2ID, "data", 4);
    waveHeader.m_subChunk2Size = dataSize;
 
    //write the header
    fwrite(&waveHeader, sizeof(SMinimalWaveFileHeader), 1, File);
 
    //write the wave data itself
    fwrite(&data[0], dataSize, 1, File);
 
    //close the file and return success
    fclose(File);
    return true;
}

bool ReadFileIntoMemory (const char *fileName, std::vector<unsigned char>& data)
{
    //open the file if we can
    FILE *file = nullptr;
    fopen_s(&file, fileName, "rb");
    if (!file)
    {
        return false;
    }

    // get the file size and resize the vector to hold the data
    fseek(file, 0, SEEK_END);
    data.resize(ftell(file));

    // read the file into the vector
    fseek(file, 0, SEEK_SET);
    fread(&data[0], 1, data.size(), file);

    // return success
    fclose(file);
    return true;
}

bool ReadWaveFile(const char *fileName, std::vector<float>& data, uint16& numChannels, uint32& sampleRate)
{
    // read the whole file into memory if we can
    std::vector<unsigned char> fileData;
    if (!ReadFileIntoMemory(fileName, fileData))
        return false;
    size_t fileIndex = 0;

	//make sure the main chunk ID is "RIFF"
	if((fileData.size() < fileIndex + 4) || memcmp(&fileData[fileIndex],"RIFF", 4))
		return false;
    fileIndex += 4;

	//get the main chunk size
	uint32 chunkSize;
    if (fileData.size() < fileIndex + 4)
        return false;
    chunkSize = *(uint32*)&fileData[fileIndex];
    fileIndex += 4;

	//make sure the format is "WAVE"
    if ((fileData.size() < fileIndex + 4) || memcmp(&fileData[fileIndex], "WAVE", 4))
        return false;
    fileIndex += 4;

    // TODO: clean up types used.

	long chunkPosFmt = -1;
	long chunkPosData = -1;
	while(chunkPosFmt == -1 || chunkPosData == -1)
	{
        // get a chunk id and chunk size if we can
        if (fileData.size() < fileIndex + 8)
            return false;

        // get the chunk id if we can
        const unsigned char* chunkID = (unsigned char*)&fileData[fileIndex];
        fileIndex += 4;
        chunkSize = *(uint32*)&fileData[fileIndex];
        fileIndex += 4;

		//if we hit a fmt
		if(!memcmp(chunkID,"fmt ", 4))
		{
			chunkPosFmt = (long)(fileIndex - 8);
		}
		//else if we hit a data
		else if(!memcmp(chunkID,"data", 4))
		{
			chunkPosData = (long)(fileIndex - 8);
		}

		//skip to the next chunk
        fileIndex += chunkSize;
	}

	//we'll use this handy struct to load in 
	SMinimalWaveFileHeader waveData;

	//load the fmt part if we can
    fileIndex = chunkPosFmt;
    if (fileData.size() < fileIndex + 24)
        return false;
    memcpy(&waveData.m_subChunk1ID, &fileData[fileIndex], 24);
    fileIndex += 24;

	//load the data part if we can
    fileIndex = chunkPosData;
    if (fileData.size() < fileIndex + 8)
        return false;
    memcpy(&waveData.m_subChunk2ID, &fileData[fileIndex], 8);
    fileIndex += 8;

	//verify a couple things about the file data
	if(waveData.m_audioFormat != 1 ||       //only pcm data
	   waveData.m_numChannels < 1 ||        //must have a channel
	   waveData.m_numChannels > 2 ||        //must not have more than 2
	   waveData.m_bitsPerSample > 32 ||     //32 bits per sample max
	   waveData.m_bitsPerSample % 8 != 0 || //must be a multiple of 8 bites
	   waveData.m_blockAlign > 8)           //blocks must be 8 bytes or lower
	{
		return false;
	}

	//figure out how many samples and blocks there are total in the source data
    int bytesPerSample = waveData.m_blockAlign / waveData.m_numChannels;
    size_t numSourceSamples = waveData.m_subChunk2Size / bytesPerSample;

	//allocate space for the source samples
    data.resize(numSourceSamples);

	//read in the source samples at whatever sample rate / number of channels it might be in
    if (fileData.size() < fileIndex + numSourceSamples * bytesPerSample)
        return false;

	for(int nIndex = 0; nIndex < numSourceSamples; ++nIndex)
	{	
        PCMToFloat(data[nIndex], &fileData[fileIndex], bytesPerSample);
        fileIndex += bytesPerSample;
	}

	//return our data
    numChannels = waveData.m_numChannels;
    sampleRate = waveData.m_sampleRate;
	return true;
}

void ChangePlaybackSpeed (const std::vector<float>& input, std::vector<float>& output, uint16 numChannels, float speedMultiplier)
{
    size_t numSrcSamples = input.size() / numChannels;
    size_t numOutSamples = (size_t)(float(numSrcSamples) / speedMultiplier);
    output.resize(numOutSamples * numChannels);

    for (size_t outSample = 0; outSample < numOutSamples; ++outSample)
    {
        float percent = float(outSample) / float(numOutSamples-1);

        float srcSampleFloat = float(numSrcSamples) * percent;
        
        size_t srcSample = size_t(srcSampleFloat);
        float srcSampleFraction = srcSampleFloat - std::floorf(srcSampleFloat);

        for (uint16 channel = 0; channel < numChannels; ++channel)
        {
            // linear interpolate samples. Cubic hermite would give better results, but this isn't the main point of the program and linear works well enough.
            float value1 = srcSample * numChannels + channel < input.size() ? input[srcSample * numChannels + channel] : *input.rbegin();
            float value2 = (srcSample + 1) * numChannels + channel < input.size() ? input[(srcSample + 1) * numChannels + channel] : *input.rbegin();
            output[outSample*numChannels + channel] = value1 * (1.0f - srcSampleFraction) + value2 * srcSampleFraction;
        }
    }
}

// SmoothStep
// a function that interpolates between 0 and 1 in a way such that it's smoother than lerp, and has flat derivatives at 0 and 1.
// https://en.wikipedia.org/wiki/Smoothstep
float SmoothStep(float x)
{
    if (x <= 0.0f)
        return 0.0f;
    else if (x >= 1.0f)
        return 1.0f;
    else
        return 3.0f * x * x - 2.0f * x * x * x;
}

// Generates a simple trapezoid type envelope (Attack, Sustain, Decay)
// Then applies smoothstep to make it more smooth
float EnvelopeGenerator (float percent, float attackDecay)
{
    if (percent < attackDecay)
        return SmoothStep(percent / attackDecay);
    else if ((1.0f - percent) < attackDecay)
        return SmoothStep((1.0f - percent) / attackDecay);
    else
        return 1.0f;
}

void GranularTimeAdjust (const std::vector<float>& input, std::vector<float>& output, uint16 numChannels, uint32 sampleRate, float timeMultiplier, float grainSizeSeconds, float envelopeSizePercent)
{
    size_t numSrcSamples = input.size() / numChannels;
    size_t numOutSamples = (size_t)(float(numSrcSamples) / timeMultiplier);
    output.resize(numOutSamples * numChannels, 0.0f);

    size_t grainSizeSamples = size_t(float(sampleRate)*grainSizeSeconds);

    size_t numGrains = numSrcSamples / grainSizeSamples;
    if (numSrcSamples % grainSizeSamples)
        numGrains++;

    size_t outSampleIndex = 0;
    for (size_t grain = 0; grain < numGrains; ++grain)
    {
        size_t grainStart = grain * grainSizeSamples;
        size_t grainEnd = (grain + 1)*grainSizeSamples;

        size_t outGrainEnd = size_t(float(grainEnd) / timeMultiplier);

        while (outSampleIndex < outGrainEnd && outSampleIndex < numOutSamples)
        {
            for (size_t sample = 0; sample < grainSizeSamples && outSampleIndex + sample < numOutSamples; ++sample)
            {
                float envelope = EnvelopeGenerator(float(sample) / float(grainSizeSamples - 1), envelopeSizePercent);
                for (size_t channel = 0; channel < numChannels; ++channel)
                {
                    output[(outSampleIndex + sample)*numChannels + channel] = input[(grainStart + sample)*numChannels + channel] * envelope;
                }
            }

            outSampleIndex += grainSizeSamples;
        }
    }
}
 
//the entry point of our application
int main(int argc, char **argv)
{
    uint16 numChannels;
    uint32 sampleRate;
    std::vector<float> source, out;
    ReadWaveFile("legend1.wav", source, numChannels, sampleRate);

#if 1
    // speed up the audio and increase pitch
    ChangePlaybackSpeed(source, out, numChannels, 1.3f);
    WriteWaveFile("out_A_Fast.wav", out, numChannels, sampleRate);

    // slow down the audio and decrease pitch
    ChangePlaybackSpeed(source, out, numChannels, 0.7f);
    WriteWaveFile("out_A_Slow.wav", out, numChannels, sampleRate);
#endif

    // speed up audio without affecting pitch
    GranularTimeAdjust(source, out, numChannels, sampleRate, 1.3f, 0.02f, 0.1f );
    WriteWaveFile("out_B_Short.wav", out, numChannels, sampleRate);

    // slow down audio without affecting pitch
    GranularTimeAdjust(source, out, numChannels, sampleRate, 0.7f, 0.01f, 0.2f );
    WriteWaveFile("out_B_Long.wav", out, numChannels, sampleRate);

    // TODO: why do we divide by the multipliers, in both GranularTimeAdjust() and ChangePlaybackSpeed()?

    // TODO: when making it shorter, do we skip granules?

    // TODO: change pitch without affecting length

    // TODO: autotune it to twinkle twinkle, and/or put it on a sine wave!
    
    // TODO: I think there may be some math issues with GranularTimeAdjust() and index calculations. I got a crash when switching to divide or passing reciprocal

    WriteWaveFile("out.wav", source, numChannels, sampleRate);
}

/*

TODO:

* this stuff needs cleaning up big time :P

* fread the entire file in at once?

* test 1, 2, 3, 4, byte formats. test their round trips too!

* actually, maybe start with source code for "lament of tim curry"
 * https://blog.demofox.org/2012/06/18/diy-synth-3-sampling-mixing-and-band-limited-wave-forms/

 * build w32 and x64

* use granular synthesis to stretch and squish a sound without affecting frequency
* also, adjust frequency without affecting length
* auto tune? with note sliding

* compare using zero crossings, cubic interpolation, and envelopes

* experiment with grain size

* make parames be apples to apples. Like... 0.7 for time adjust is 1.3 for playback speed adjust?? make it take the same values

BLOG:

* Explain algorithm, give results, link to code and also include it!

https://granularsynthesis.com/guide.php

Read these:
https://www.soundonsound.com/techniques/granular-synthesis
https://www.granularsynthesis.com/hthesis/grain.html

* note the thing about how on graphics cards there are two -1 values, to make conversion between types easier. same issues as pcm <-> float here!
 * http://www.yosoygames.com.ar/wp/2018/03/vertex-formats-part-1-compression/

*/