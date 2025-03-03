/*
Description: 解析read的name中的信息，比如tile, x, y等

Copyright : All right reserved by ICT

Author : Zhang Zhonghai
Date : 2023/11/6
*/
#ifndef READ_NAME_PARSER_H_
#define READ_NAME_PARSER_H_

#include <spdlog/spdlog.h>
#include <stdint.h>

#include <regex>
#include <string>

#include "read_ends.h"

// using std::regex;
using std::cmatch;
using std::regex;
using std::string;



/**
 * Provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile
 * should only allow non-zero positive integers, x and y coordinates may be
 * negative. 非线程安全
 */
struct ReadNameParser {
    /**
     * The read name regular expression (regex) is used to extract three pieces
     * of information from the read name: tile, x location, and y location.  Any
     * read name regex should parse the read name to produce these and only
     * these values.  An example regex is:
     *  (?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$
     * which assumes that fields in the read name are delimited by ':' and the
     * last three fields correspond to the tile, x and y locations, ignoring any
     * trailing non-digit characters.
     *
     * The default regex is optimized for fast parsing (see {@link
     * #getLastThreeFields(String, char, int[])}) by searching for the last
     * three fields, ignoring any trailing non-digit characters, assuming the
     * delimiter ':'.  This should consider correctly read names where we have 5
     * or 7 field with the last three fields being tile/x/y, as is the case for
     * the majority of read names produced by Illumina technology.
     */
    const string DEFAULT_READ_NAME_REGEX = "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$";
    bool warnAboutRegexNotMatching = true;

    static bool sWrongNameFormat;

    string readNameStored = "";
    PhysicalLocation physicalLocationStored;
    int tmpLocationFields[3];                // for optimization of addLocationInformation
    bool useOptimizedDefaultParsing = true;  // was the regex default?
    string readNameRegex = DEFAULT_READ_NAME_REGEX;
    regex readNamePattern;

    ReadNameParser() : ReadNameParser(DEFAULT_READ_NAME_REGEX) {}
    ReadNameParser(const string &strReadNameRegex) : ReadNameParser(strReadNameRegex, true) {}
    ReadNameParser(const string &strReadNameRegex, bool isWarn) {
        readNameRegex = strReadNameRegex;
        if (strReadNameRegex == DEFAULT_READ_NAME_REGEX)
            useOptimizedDefaultParsing = true;
        else
            useOptimizedDefaultParsing = false;
        readNamePattern = std::regex(strReadNameRegex, std::regex_constants::optimize);
        warnAboutRegexNotMatching = isWarn;
    }

    /* 重新设置readNameRegex */
    void SetReadNameRegex(const string &strReadNameRegex) {
        readNameRegex = strReadNameRegex;
        if (strReadNameRegex == DEFAULT_READ_NAME_REGEX)
            useOptimizedDefaultParsing = true;
        else
            useOptimizedDefaultParsing = false;
        readNamePattern = std::regex(strReadNameRegex, std::regex_constants::optimize);
        // readNamePattern = strReadNameRegex;
    }

    /* 添加测序时候的tile x y 信息 */
    bool AddLocationInformation(const string &readName, PhysicalLocation *loc) {
        if (!(readName == readNameStored)) {
            if (ReadLocationInformation(readName, loc)) {
                readNameStored = readName;
                physicalLocationStored = *loc;
                return true;
            }
            // return false if read name cannot be parsed
            return false;
        } else {
            *loc = physicalLocationStored;
            return true;
        }
    }

    /**
     * Method used to extract tile/x/y from the read name and add it to the
     * PhysicalLocationShort so that it can be used later to determine optical
     * duplication
     *
     * @param readName the name of the read/cluster
     * @param loc the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form,
     * false otherwise
     */
    bool ReadLocationInformation(const string &readName, PhysicalLocation *loc) {
        try {
            // Optimized version if using the default read name regex (== used on purpose):
            if (useOptimizedDefaultParsing) {
                const int fields = getLastThreeFields(readName, ':');
                if (!(fields == 5 || fields == 7)) {
                    if (warnAboutRegexNotMatching) {
                        spdlog::warn(
                            "Default READ_NAME_REGEX '{}' did not match read "
                            "name '{}'."
                            "You may need to specify a READ_NAME_REGEX in "
                            "order to correctly identify optical duplicates.  "
                            "Note that this message will not be emitted again "
                            "even if other read names do not match the regex.",
                            readNameRegex.c_str(), readName.c_str());
                        warnAboutRegexNotMatching = false;
                        sWrongNameFormat = true;
                    }
                    return false;
                }
                loc->tile = (int16_t)tmpLocationFields[0];
                loc->x = tmpLocationFields[1];
                loc->y = tmpLocationFields[2];
                return true;
            } else if (readNameRegex.empty()) {
                return false;
            } else {
                // Standard version that will use the regex
                cmatch m;
                if (std::regex_match(readName.c_str(), m, readNamePattern)) {
                    loc->tile = std::stoi(m[1].str());
                    loc->x = std::stoi(m[2].str());
                    loc->y = std::stoi(m[3].str());
                    return true;
                } else {
                    if (warnAboutRegexNotMatching) {
                        spdlog::warn(
                            "READ_NAME_REGEX '{}' did not match read name '{}'."
                            "Your regex may not be correct.  "
                            "Note that this message will not be emitted again "
                            "even if other read names do not match the regex.",
                            readNameRegex.c_str(), readName.c_str());
                        warnAboutRegexNotMatching = false;
                        sWrongNameFormat = true;
                    }
                    return false;
                }
            }
        } catch (const std::runtime_error &e) {
            if (warnAboutRegexNotMatching) {
                spdlog::warn(
                    "A field parsed out of a read name was expected to contain "
                    "an integer and did not. READ_NAME_REGEX: {}; Read name: "
                    "{}; Error Msg: {}",
                    readNameRegex.c_str(), readName.c_str(), e.what());
                warnAboutRegexNotMatching = false;
                sWrongNameFormat = true;
            }
        } catch (...) {
            if (warnAboutRegexNotMatching) {
                spdlog::warn(
                    "A field parsed out of a read name was expected to contain "
                    "an integer and did not. READ_NAME_REGEX: {}; Read name: "
                    "{}",
                    readNameRegex.c_str(), readName.c_str());
                warnAboutRegexNotMatching = false;
                sWrongNameFormat = true;
            }
        }

        return true;
    }

    /**
     * Given a string, splits the string by the delimiter, and returns the the
     * last three fields parsed as integers.  Parsing a field considers only a
     * sequence of digits up until the first non-digit character.  The three
     * values are stored in the passed-in array.
     *
     * @throws NumberFormatException if any of the tokens that should contain
     * numbers do not start with parsable numbers
     */
    int getLastThreeFields(const string &readName, char delim) {
        int tokensIdx = 2;  // start at the last token
        int numFields = 0;
        int i, endIdx;
        endIdx = readName.size();
        // find the last three tokens only
        for (i = (int)readName.size() - 1; 0 <= i && 0 <= tokensIdx; i--) {
            if (readName.at(i) == delim || 0 == i) {
                numFields++;
                const int startIdx = (0 == i) ? 0 : (i + 1);
                tmpLocationFields[tokensIdx] = std::stoi(readName.substr(startIdx, endIdx - startIdx));
                tokensIdx--;
                endIdx = i;
            }
        }
        // continue to find the # of fields
        while (0 <= i) {
            if (readName.at(i) == delim || 0 == i)
                numFields++;
            i--;
        }
        if (numFields < 3) {
            tmpLocationFields[0] = tmpLocationFields[1] = tmpLocationFields[2] = -1;
            numFields = -1;
        }

        return numFields;
    }
};

#endif