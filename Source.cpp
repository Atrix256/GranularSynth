#include <stdio.h>
#include <memory.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
 
// typedefs
typedef uint16_t    uint16;
typedef uint32_t    uint32;
typedef int32_t     int32;

const float c_pi = 3.14159265359f;
 
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

enum class ECrossFade
{
    None,
    In,
    Out,
};

inline void FloatToPCM(unsigned char *PCM, const float& in, size_t numBytes)
{
    // 8 bit is unsigned
    if (numBytes == 1)
    {
        PCM[0] = unsigned char((in * 0.5f + 0.5f) * 255.0f);
        return;
    }

    // casting to double because floats can't exactly store 0x7fffffff, but doubles can.
    // Details of that: https://blog.demofox.org/2017/11/21/floating-point-precision/
    uint32 data;
    if (in < 0.0f)
        data = uint32(double(in) * double(0x80000000));
    else
        data = uint32(double(in) * double(0x7fffffff));

    switch (numBytes)
    {
        case 4: PCM[3] = ((data >> 24) & 0xFF); PCM[2] = ((data >> 16) & 0xFF); PCM[1] = ((data >> 8) & 0xFF); PCM[0] = (data & 0xFF); break;
        case 3: PCM[2] = ((data >> 24) & 0xFF); PCM[1] = ((data >> 16) & 0xFF); PCM[0] = ((data >> 8) & 0xFF); break;
        case 2: PCM[1] = ((data >> 24) & 0xFF); PCM[0] = ((data >> 16) & 0xFF); break;
    }
}

inline void PCMToFloat(float& out, const unsigned char *PCM, size_t numBytes)
{
    // 8 bit is unsigned
    if (numBytes == 1)
    {
        out = (float(PCM[0]) / float(255.0f)) * 2.0f - 1.0f;
        return;
    }

    uint32 data = 0;
    switch (numBytes)
    {
        case 4: data = (uint32(PCM[3]) << 24) | (uint32(PCM[2]) << 16) | (uint32(PCM[1]) << 8) | uint32(PCM[0]); break;
        case 3: data = (uint32(PCM[2]) << 24) | (uint32(PCM[1]) << 16) | (uint32(PCM[0]) << 8); break;
        case 2: data = (uint32(PCM[1]) << 24) | (uint32(PCM[0]) << 16); break;
    }

    // casting to double because floats can't exactly store 0x7fffffff, but doubles can.
    // Details of that: https://blog.demofox.org/2017/11/21/floating-point-precision/
    if (data & 0x80000000)
        out = float(double(int32(data)) / double(0x80000000));
    else
        out = float(double(data) / double(0x7fffffff));
}

// numBytes can be 1, 2, 3, or 4.
// Coresponding to 8 bit, 16 bit, 24 bit, and 32 bit audio.
bool WriteWaveFile(const char *fileName, std::vector<float>& dataFloat, uint16 numChannels, uint32 sampleRate, uint16 numBytes)
{
    std::vector<unsigned char> data;
    data.resize(dataFloat.size() * numBytes);
    for (size_t i = 0; i < dataFloat.size(); ++i)
        FloatToPCM((unsigned char*)&data[i*numBytes], dataFloat[i], numBytes);

    uint32 dataSize = (uint32)data.size();
    uint16 bitsPerSample = numBytes * 8;
 
    //open the file if we can
    FILE *File = nullptr;
    fopen_s(&File, fileName, "w+b");
    if (!File)
    {
        printf("[-----ERROR-----] Could not open %s for writing.\n", fileName);
        return false;
    }
 
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
    printf("%s saved.\n", fileName);
    return true;
}

bool ReadFileIntoMemory (const char *fileName, std::vector<unsigned char>& data)
{
    //open the file if we can
    FILE *file = nullptr;
    fopen_s(&file, fileName, "rb");
    if (!file)
    {
        printf("[-----ERROR-----]Could not open %s for reading.\n", fileName);
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

bool ReadWaveFile(const char *fileName, std::vector<float>& data, uint16& numChannels, uint32& sampleRate, uint16& numBytes)
{
    // read the whole file into memory if we can
    std::vector<unsigned char> fileData;
    if (!ReadFileIntoMemory(fileName, fileData))
        return false;
    size_t fileIndex = 0;

	//make sure the main chunk ID is "RIFF"
    if ((fileData.size() < fileIndex + 4) || memcmp(&fileData[fileIndex], "RIFF", 4))
    {
        printf("[-----ERROR-----]%s is an invalid input file. (1)\n", fileName);
        return false;
    }
    fileIndex += 4;

	//get the main chunk size
	uint32 chunkSize;
    if (fileData.size() < fileIndex + 4)
    {
        printf("[-----ERROR-----]%s is an invalid input file. (2)\n", fileName);
        return false;
    }
    chunkSize = *(uint32*)&fileData[fileIndex];
    fileIndex += 4;

	//make sure the format is "WAVE"
    if ((fileData.size() < fileIndex + 4) || memcmp(&fileData[fileIndex], "WAVE", 4))
    {
        printf("[-----ERROR-----]%s is an invalid input file. (3)\n", fileName);
        return false;
    }
    fileIndex += 4;

	size_t chunkPosFmt = -1;
	size_t chunkPosData = -1;
	while(chunkPosFmt == -1 || chunkPosData == -1)
	{
        // get a chunk id and chunk size if we can
        if (fileData.size() < fileIndex + 8)
        {
            printf("[-----ERROR-----]%s is an invalid input file. (4)\n", fileName);
            return false;
        }

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
    {
        printf("[-----ERROR-----]%s is an invalid input file. (5)\n", fileName);
        return false;
    }
    memcpy(&waveData.m_subChunk1ID, &fileData[fileIndex], 24);
    fileIndex += 24;

	//load the data part if we can
    fileIndex = chunkPosData;
    if (fileData.size() < fileIndex + 8)
    {
        printf("[-----ERROR-----]%s is an invalid input file. (6)\n", fileName);
        return false;
    }
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
        printf("[-----ERROR-----]%s is an invalid input file. (7)\n", fileName);
        return false;
    }

	//figure out how many samples and blocks there are total in the source data
    size_t bytesPerSample = waveData.m_blockAlign / waveData.m_numChannels;
    size_t numSourceSamples = waveData.m_subChunk2Size / bytesPerSample;

	//allocate space for the source samples
    data.resize(numSourceSamples);

	//read in the source samples at whatever sample rate / number of channels it might be in
    if (fileData.size() < fileIndex + numSourceSamples * bytesPerSample)
    {
        printf("[-----ERROR-----]%s is an invalid input file. (8)\n", fileName);
        return false;
    }

	for(size_t nIndex = 0; nIndex < numSourceSamples; ++nIndex)
	{	
        PCMToFloat(data[nIndex], &fileData[fileIndex], bytesPerSample);
        fileIndex += bytesPerSample;
	}

	//return our data
    numChannels = waveData.m_numChannels;
    sampleRate = waveData.m_sampleRate;
    numBytes = waveData.m_bitsPerSample / 8;

    printf("%s loaded.\n", fileName);
	return true;
}

inline float SampleChannelFractional (const std::vector<float>& input, float sampleFloat, uint16 channel, uint16 numChannels)
{
    // TODO: cubic hermite interpolate for better results

    size_t sample = size_t(sampleFloat);
    float sampleFraction = sampleFloat - std::floorf(sampleFloat);

    float value1 = sample * numChannels + channel < input.size() ? input[sample * numChannels + channel] : *input.rbegin();
    float value2 = (sample + 1) * numChannels + channel < input.size() ? input[(sample + 1) * numChannels + channel] : *input.rbegin();
    return value1 * (1.0f - sampleFraction) + value2 * sampleFraction;
}

// Resample
void TimeAdjust (const std::vector<float>& input, std::vector<float>& output, uint16 numChannels, float timeMultiplier)
{
    size_t numSrcSamples = input.size() / numChannels;
    size_t numOutSamples = (size_t)(float(numSrcSamples) * timeMultiplier);
    output.resize(numOutSamples * numChannels);

    for (size_t outSample = 0; outSample < numOutSamples; ++outSample)
    {
        float percent = float(outSample) / float(numOutSamples-1);

        float srcSampleFloat = float(numSrcSamples) * percent;

        for (uint16 channel = 0; channel < numChannels; ++channel)
            output[outSample*numChannels + channel] = SampleChannelFractional(input, srcSampleFloat, channel, numChannels);
    }
}

// writes a grain to the output buffer, applying a fade in or fade out at the beginning if it should.
size_t SplatGrainToOutput (const std::vector<float>& input, std::vector<float>& output, uint16 numChannels, size_t grainStart, size_t grainSize, size_t outputSampleIndex, ECrossFade crossFade, size_t crossFadeSize, float pitchMultiplier)
{
    // Enforce input and output array size constraints
    size_t numInputSamples = input.size() / numChannels;
    if (grainStart + grainSize > numInputSamples)
        grainSize = numInputSamples - grainStart;

    size_t numOutputSamples = output.size() / numChannels;
    if (outputSampleIndex + grainSize > numOutputSamples)
        grainSize = numOutputSamples - outputSampleIndex;

    // calculate starting indices
    size_t inputIndex = grainStart * numChannels;
    size_t outputIndex = outputSampleIndex * numChannels;

    // write the samples
    for (size_t sample = 0; sample < grainSize; ++sample)
    {
        // calculate envelope for this sample
        float envelope = 1.0f;
        if (crossFade != ECrossFade::None)
        {
            if (sample <= crossFadeSize)
                envelope = float(sample) / float(crossFadeSize);
            if (crossFade == ECrossFade::Out)
                envelope = 1.0f - envelope;
        }

        // write the enveloped sample
        for (size_t channel = 0; channel < numChannels; ++channel)
            output[outputIndex + channel] += input[inputIndex + channel] * envelope;

        // move to the next samples
        inputIndex += numChannels;
        outputIndex += numChannels;
    }

    // return how many samples we wrote
    return grainSize;
}

void GranularTimeAdjust (const std::vector<float>& input, std::vector<float>& output, uint16 numChannels, uint32 sampleRate, float timeMultiplier, float pitchMultiplier, float grainSizeSeconds, float crossFadeSeconds)
{
    // calculate size of output buffer and resize it
    size_t numInputSamples = input.size() / numChannels;
    size_t numOutputSamples = (size_t)(float(numInputSamples) * timeMultiplier);
    output.clear();
    output.resize(numOutputSamples * numChannels, 0.0f);

    // calculate how many grains are in the input data
    size_t grainSizeSamples = size_t(float(sampleRate)*grainSizeSeconds);
    size_t numGrains = numInputSamples / grainSizeSamples;
    if (numInputSamples % grainSizeSamples)
        numGrains++;

    // calculate the cross fade size
    size_t crossFadeSizeSamples = size_t(float(sampleRate)*crossFadeSeconds);

    // Repeat each grain 0 or more times to make the output be the correct size
    size_t outputSampleIndex = 0;
    size_t lastGrainWritten = -1;
    for (size_t grain = 0; grain < numGrains; ++grain)
    {
        // calculate the boundaries of the grain
        size_t inputGrainStart = grain * grainSizeSamples;

        // calculate the end of where this grain should go in the output buffer
        size_t outputSampleWindowEnd = size_t(float(inputGrainStart + grainSizeSamples) * timeMultiplier);

        // Splat out zero or more copies of the grain to get our output to be at least as far as we want it to be.
        // Zero copies happens when we shorten time and need to cut pieces (grains) out of the original sound
        while (outputSampleIndex < outputSampleWindowEnd)
        {
            // if we are writing our first grain, or the last grain we wrote was the previous grain, then we don't need to do a cross fade`
            if ((lastGrainWritten == -1) || (lastGrainWritten == grain - 1))
            {
                outputSampleIndex += SplatGrainToOutput(input, output, numChannels, inputGrainStart, grainSizeSamples, outputSampleIndex, ECrossFade::None, crossFadeSizeSamples, pitchMultiplier);
                lastGrainWritten = grain;
                continue;
            }

            // else we need to fade out the old grain and then fade in the new one.
            // NOTE: fading out the old grain means starting to play the grain after the last one and bringing it's volume down to zero.
            SplatGrainToOutput(input, output, numChannels, (lastGrainWritten + 1) * grainSizeSamples, grainSizeSamples, outputSampleIndex, ECrossFade::Out, crossFadeSizeSamples, pitchMultiplier);
            outputSampleIndex += SplatGrainToOutput(input, output, numChannels, inputGrainStart, grainSizeSamples, outputSampleIndex, ECrossFade::In, crossFadeSizeSamples, pitchMultiplier);
            lastGrainWritten = grain;
        }
    }
}

template <typename LAMBDA>
void GranularTimeAdjustDynamic (const std::vector<float>& input, std::vector<float>& output, uint16 numChannels, uint32 sampleRate, float grainSizeSeconds, float crossFadeSeconds, const LAMBDA& timeMultiplierCallback )
{
    // calculate how many grains are in the input data
    size_t numInputSamples = input.size() / numChannels;
    size_t grainSizeSamples = size_t(float(sampleRate)*grainSizeSeconds);
    size_t numGrains = numInputSamples / grainSizeSamples;
    if (numInputSamples % grainSizeSamples)
        numGrains++;

    // calculate size of output buffer and resize it
    size_t numOutputSamples = 0;
    for (size_t i = 0; i < numGrains; ++i)
    {
        size_t grainStart = i * grainSizeSamples;
        size_t grainEnd = grainStart + grainSizeSamples;
        grainEnd = std::min(grainEnd, input.size());
        size_t grainSize = grainEnd - grainStart;

        float percent = float(i) / float(numGrains);
        float timeMultiplier = timeMultiplierCallback(percent);

        numOutputSamples += (size_t)(float(grainSize) * timeMultiplier);
    }
    output.clear();
    output.resize(numOutputSamples * numChannels, 0.0f);

    // calculate the cross fade size
    size_t crossFadeSizeSamples = size_t(float(sampleRate)*crossFadeSeconds);

    // Repeat each grain 0 or more times to make the output be the correct size
    size_t outputSampleIndex = 0;
    size_t lastGrainWritten = -1;
    size_t outputSampleWindowEnd = 0;
    for (size_t grain = 0; grain < numGrains; ++grain)
    {
        // calculate the boundaries of the grain
        size_t inputGrainStart = grain * grainSizeSamples;

        // calculate the end of where this grain should go in the output buffer
        float percent = float(grain) / float(numGrains);
        float timeMultiplier = timeMultiplierCallback(percent);
        outputSampleWindowEnd += size_t(float(grainSizeSamples) * timeMultiplier);

        // Splat out zero or more copies of the grain to get our output to be at least as far as we want it to be.
        // Zero copies happens when we shorten time and need to cut pieces (grains) out of the original sound
        while (outputSampleIndex < outputSampleWindowEnd)
        {
            // if we are writing our first grain, or the last grain we wrote was the previous grain, then we don't need to do a cross fade`
            if ((lastGrainWritten == -1) || (lastGrainWritten == grain - 1))
            {
                outputSampleIndex += SplatGrainToOutput(input, output, numChannels, inputGrainStart, grainSizeSamples, outputSampleIndex, ECrossFade::None, crossFadeSizeSamples);
                lastGrainWritten = grain;
                continue;
            }

            // else we need to fade out the old grain and then fade in the new one.
            // NOTE: fading out the old grain means starting to play the grain after the last one and bringing it's volume down to zero.
            SplatGrainToOutput(input, output, numChannels, (lastGrainWritten + 1) * grainSizeSamples, grainSizeSamples, outputSampleIndex, ECrossFade::Out, crossFadeSizeSamples);
            outputSampleIndex += SplatGrainToOutput(input, output, numChannels, inputGrainStart, grainSizeSamples, outputSampleIndex, ECrossFade::In, crossFadeSizeSamples);
            lastGrainWritten = grain;
        }
    }
}

//the entry point of our application
int main(int argc, char **argv)
{
    // load the wave file
    uint16 numChannels;
    uint32 sampleRate;
    uint16 numBytes;
    std::vector<float> source, out, out2, sourceLeft, sourceRight;
    ReadWaveFile("legend1.wav", source, numChannels, sampleRate, numBytes);

    // split the stereo into two mono channels
    sourceLeft.resize(source.size() / numChannels);
    sourceRight.resize(sourceLeft.size());
    for (size_t i = 0; i < sourceLeft.size(); ++i)
    {
        sourceLeft[i] = source[i*numChannels];
        sourceRight[i] = source[i*numChannels + 1];
    }

#if 1
    // speed up the audio and increase pitch
    {
        TimeAdjust(source, out, numChannels, 0.7f);
        WriteWaveFile("out_A_FastHigh.wav", out, numChannels, sampleRate, numBytes);

        TimeAdjust(source, out, numChannels, 0.4f);
        WriteWaveFile("out_A_FasterHigher.wav", out, numChannels, sampleRate, numBytes);
    }

    // slow down the audio and decrease pitch
    {
        TimeAdjust(source, out, numChannels, 1.3f);
        WriteWaveFile("out_A_SlowLow.wav", out, numChannels, sampleRate, numBytes);

        TimeAdjust(source, out, numChannels, 2.1f);
        WriteWaveFile("out_A_SlowerLower.wav", out, numChannels, sampleRate, numBytes);
    }

    // speed up audio without affecting pitch
    {
        GranularTimeAdjust(source, out, numChannels, sampleRate, 0.7f, 1.0f, 0.02f, 0.002f);
        WriteWaveFile("out_B_Fast.wav", out, numChannels, sampleRate, numBytes);

        GranularTimeAdjust(source, out, numChannels, sampleRate, 0.4f, 1.0f, 0.02f, 0.002f);
        WriteWaveFile("out_B_Faster.wav", out, numChannels, sampleRate, numBytes);
    }

    // slow down audio without affecting pitch
    {
        GranularTimeAdjust(source, out, numChannels, sampleRate, 1.3f, 1.0f, 0.02f, 0.002f);
        WriteWaveFile("out_B_Slow.wav", out, numChannels, sampleRate, numBytes);

        GranularTimeAdjust(source, out, numChannels, sampleRate, 2.1f, 1.0f, 0.02f, 0.002f);
        WriteWaveFile("out_B_Slower.wav", out, numChannels, sampleRate, numBytes);
    }

    // Make pitch higher without affecting length
    {
        GranularTimeAdjust(source, out2, numChannels, sampleRate, 1.0f / 0.7f, 1.0f, 0.02f, 0.002f);
        TimeAdjust(out2, out, numChannels, 0.7f);
        WriteWaveFile("out_C_High.wav", out, numChannels, sampleRate, numBytes);

        GranularTimeAdjust(source, out2, numChannels, sampleRate, 1.0f / 0.4f, 1.0f, 0.02f, 0.002f);
        TimeAdjust(out2, out, numChannels, 0.4f);
        WriteWaveFile("out_C_Higher.wav", out, numChannels, sampleRate, numBytes);
    }

    // make pitch lower without affecting length
    {
        GranularTimeAdjust(source, out2, numChannels, sampleRate, 1.0f / 1.3f, 1.0f, 0.02f, 0.002f);
        TimeAdjust(out2, out, numChannels, 1.3f);
        WriteWaveFile("out_C_Low.wav", out, numChannels, sampleRate, numBytes);

        GranularTimeAdjust(source, out2, numChannels, sampleRate, 1.0f / 2.1f, 1.0f, 0.02f, 0.002f);
        TimeAdjust(out2, out, numChannels, 2.1f);
        WriteWaveFile("out_C_Lower.wav", out, numChannels, sampleRate, numBytes);
    }
#endif

    /*
    // TODO: this more adjusts playback speed on a sine wave, not pitch!!  I think you would need a TimeAdjustDynamic() to make it work
    // TODO: may need some kind of normalization constant (calculate it!) to get this to integrate to 1.0 to not affect pitch or speed (whichever should be left alone)

    // TODO: i think the function has to integrate to 1.
    // Note that this does. Check wolfram alpha: integrate y=(sin(2*pi*x)+1.0) from 0 to 1

    // adjust pitch on a sine wave
    GranularTimeAdjustDynamic(source, out2, numChannels, sampleRate, 0.02f, 0.002f,
        [] (float percent)
        {
            // TODO: i think whatever function you give it needs to integrate to 1.0 from 0 to 1. maybe get a normalization constant when calculating length!
            return std::sinf(sin(2 * c_pi*percent) + 1.0);

            //return (std::sinf(percent * c_pi * 10.0f) * 0.5f + 0.5f) * 1.5f + 0.5f;
        }
    );
    float ratio = float(source.size()) / float(out2.size());
    TimeAdjust(out2, out, numChannels, ratio);
    WriteWaveFile("out_D.wav", out, numChannels, sampleRate, numBytes);
    */

    system("pause");

    // TODO: do one where the speed is slow but the pitch is high, and one where the speed is fast but the pitch is low
}

/*

TODO:

* do a zero crossings version
* for zero crossings one, do each channel individually yeah! split the stereo into two mono's, process them, and recombine!

* get cubic hermite interpolation working instead of lerp for sample interpolation, for better results

* experiment with grain size and envelope size etc

BLOG:

* Explain algorithm, give results, link to code and also include it!
* get links & info from email

https://granularsynthesis.com/guide.php

Read these:
https://www.soundonsound.com/techniques/granular-synthesis
https://www.granularsynthesis.com/hthesis/grain.html

* note the thing about how on graphics cards there are two -1 values, to make conversion between types easier. same issues as pcm <-> float here!
 * http://www.yosoygames.com.ar/wp/2018/03/vertex-formats-part-1-compression/

* zero crossings don't work in stereo unless you handle each channel separately

*/