/*
Description: Murmur哈希

Copyright : All right reserved by ICT

Author : Zhang Zhonghai
Date : 2023/11/6
*/

#include <random>
#include <string>

using std::string;

/**
 * Provides an implementation of the Murmur3_32 hash algorithm that has
 * desirable properties in terms of randomness and uniformity of the
 * distribution of output values that make it a useful hashing algorithm for
 * downsampling.
 */

struct Murmur3 {
    int seed_ = 0;
    /** Hashes a character stream to an int using Murmur3. */
    int HashUnencodedChars(const string &input) {
        int h1 = this->seed_;

        // step through the CharSequence 2 chars at a time
        const int length = input.size();
        for (int i = 1; i < length; i += 2) {
            int k1 = input.at(i - 1) | (input.at(i) << 16);
            k1 = mixK1(k1);
            h1 = mixH1(h1, k1);
        }

        // deal with any remaining characters
        if ((length & 1) == 1) {
            int k1 = input.at(length - 1);
            k1 = mixK1(k1);
            h1 ^= k1;
        }

        return fmix(h1, 2 * length);
    }

    static Murmur3 &Instance() {
        static Murmur3 instance;
        return instance;
    }

    static int mixK1(int k1) {
        const int c1 = 0xcc9e2d51;
        const int c2 = 0x1b873593;
        k1 *= c1;
        k1 = k1 << 15;
        k1 *= c2;
        return k1;
    }

    static int mixH1(int h1, int k1) {
        h1 ^= k1;
        h1 = h1 << 13;
        h1 = h1 * 5 + 0xe6546b64;
        return h1;
    }

    // Finalization mix - force all bits of a hash block to avalanche
    static int fmix(int h1, int length) {
        h1 ^= length;
        h1 ^= (unsigned int)h1 >> 16;
        h1 *= 0x85ebca6b;
        h1 ^= (unsigned int)h1 >> 13;
        h1 *= 0xc2b2ae35;
        h1 ^= (unsigned int)h1 >> 16;
        return h1;
    }

private:
    Murmur3() {
        auto &&rd = std::random_device{};
        seed_ = rd();
    }
};