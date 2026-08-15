//
// Copyright 2022 Eugene Gershnik
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file or at
// https://github.com/gershnik/argum/blob/master/LICENSE
//
#ifndef HEADER_ARGUM_ARGUM_H_INCLUDED
#define HEADER_ARGUM_ARGUM_H_INCLUDED


#include <algorithm>
#include <assert.h>
#include <charconv>
#include <concepts>
#include <ctype.h>
#include <exception>
#include <filesystem>
#include <functional>
#ifdef _WIN32
    #include <io.h>
#endif
#include <iterator>
#include <limits>
#include <limits.h>
#include <map>
#include <math.h>
#include <memory>
#include <optional>
#include <regex>
#include <span>
#include <stack>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <string_view>
#if !defined(_WIN32) && __has_include(<sys/ioctl.h>)
    #include <sys/ioctl.h>
#endif
#include <system_error>
#if !defined(_WIN32) && __has_include(<termios.h>)
    #include <termios.h>
#endif
#include <tuple>
#include <type_traits>
#if !defined(_WIN32) && __has_include(<unistd.h>)
    #include <unistd.h>
#endif
#include <variant>
#include <vector>
#include <wchar.h>
#include <wctype.h>
#ifdef _WIN32
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #ifndef WIN32_LEAN_AND_MEAN
        #define WIN32_LEAN_AND_MEAN
    #endif
    #include <Windows.h>
    #ifdef min
        #undef min
    #endif
    #ifdef max
        #undef max
    #endif
#endif

#ifndef HEADER_ARGUM_UBMRELLA_INCLUDED
#define HEADER_ARGUM_UBMRELLA_INCLUDED

#ifndef HEADER_ARGUM_PARSER_H_INCLUDED
#define HEADER_ARGUM_PARSER_H_INCLUDED

#ifndef HEADER_ARGUM_CHAR_CONSTANTS_H_INCLUDED
#define HEADER_ARGUM_CHAR_CONSTANTS_H_INCLUDED

#ifndef HEADER_ARGUM_COMMON_H_INCLUDED
#define HEADER_ARGUM_COMMON_H_INCLUDED



#ifndef ARGUM_MOD_EXPORTED
    #define ARGUM_MOD_EXPORTED
#endif

#define ARGUM_DECLARE_FRIENDLY_NAME(stem, type, prefix) ARGUM_MOD_EXPORTED using prefix ## stem = Basic##stem<type>;
    
#define ARGUM_DECLARE_FRIENDLY_NAMES(stem) \
    ARGUM_DECLARE_FRIENDLY_NAME(stem, char, ) \
    ARGUM_DECLARE_FRIENDLY_NAME(stem, wchar_t, W)

#define ARGUM_CONCAT1(a, b) a ## b
#define ARGUM_CONCAT(a, b) ARGUM_CONCAT1(a, b)

#ifdef __COUNTER__
    #define ARGUM_UNIQUE_NAME(stem) ARGUM_CONCAT(stem, __COUNTER__)
#else
    #define ARGUM_UNIQUE_NAME(stem) ARGUM_CONCAT(stem, __LINE__)
#endif

#ifndef NDEBUG
    #define ARGUM_ALWAYS_ASSERT(x) assert(x)
#else
    #define ARGUM_ALWAYS_ASSERT(x)  ((void) ((x) ? ((void)0) : ::Argum::terminateApplication("failed assertion `" #x "'")))
#endif

#define ARGUM_INVALID_ARGUMENT(message) ::Argum::terminateApplication(message)


#if defined(__GNUC__)
    #ifndef __EXCEPTIONS
        #define ARGUM_NO_THROW
    #endif
#elif defined(_MSC_VER)
    #if !defined(_CPPUNWIND) || _HAS_EXCEPTIONS == 0
        #define ARGUM_NO_THROW
    #endif
#endif

#ifndef ARGUM_NO_THROW
    #define ARGUM_RAISE_EXCEPTION(x) throw x
#else
    #ifndef ARGUM_USE_EXPECTED
        #define ARGUM_USE_EXPECTED 1
    #endif
    #define ARGUM_RAISE_EXCEPTION(x) ::Argum::terminateApplication((x).what())
#endif 

namespace Argum {

    [[noreturn]] auto terminateApplication(const char * message) -> void; 

    #ifndef ARGUM_CUSTOM_TERMINATE
        [[noreturn]] inline auto terminateApplication(const char * message) -> void { 
            fprintf(stderr, "%s\n", message); 
            fflush(stderr); 
            #ifndef NDEBUG
                assert(false);
            #endif
            std::terminate();
        }
    #endif

    template<class Char>
    concept Character = std::is_same_v<Char, char> || 
                        std::is_same_v<Char, wchar_t>;

    template<class X>
    concept StringLike = requires(X name) {
        requires Character<std::remove_cvref_t<decltype(name[0])>>;
        requires std::is_convertible_v<X, std::basic_string_view<std::remove_cvref_t<decltype(name[0])>>>;
    };

    template<StringLike S> using CharTypeOf = std::remove_cvref_t<decltype(std::declval<S>()[0])>;

    template<class X, class Char>
    concept StringLikeOf = Character<Char> && std::is_convertible_v<X, std::basic_string_view<Char>>;


    template<class T, class Char>
    concept ArgIterator  = requires(T & t) {
        t + 0u;
        t - 0u;
        t++;
        t--;
        std::is_convertible_v<decltype(*t), std::basic_string_view<Char>>;
        std::is_same_v<decltype(t != t), bool>;
        std::is_same_v<decltype(t == t), bool>;
        std::is_convertible_v<decltype(t[0]), std::basic_string_view<Char>>;
    };

    template<class T, class Char>
    concept ArgRange = requires(T & t) {
        std::begin(t);
        requires ArgIterator<decltype(std::begin(t)), Char>;
        std::end(t);
        requires std::is_same_v<decltype(std::begin(t) != std::end(t)), bool>;
        requires std::is_convertible_v<decltype(t[0]), std::basic_string_view<Char>>;
    };


    template<Character Char>
    auto makeString(Char c) -> std::basic_string<Char> {
        return std::basic_string<Char>(1, c);
    }

    template<Character Char>
    auto makeString(std::basic_string_view<Char> view) -> std::basic_string<Char> {
        return std::basic_string<Char>(view);
    }

    template<Character Char>
    auto makeString(const Char * str) -> std::basic_string<Char> {
        return std::basic_string<Char>(str);
    }

    template<Character Char>
    auto makeString(std::basic_string<Char> && str) -> std::basic_string<Char> {
        return std::basic_string<Char>(std::move(str));
    }

    template<Character Char>
    auto makeString(const std::basic_string<Char> & str) -> std::basic_string<Char> {
        return std::basic_string<Char>(str);
    }

    template<class T, class Char>
    concept StringConvertibleOf = Character<Char> && requires(T && t) {
        {makeString(t)} -> std::same_as<std::basic_string<Char>>;
    };

    ARGUM_MOD_EXPORTED
    template<class It1, class It2, class Joiner>
    auto join(It1 first, It2 last, Joiner joiner) -> decltype(*first + joiner) {
        
        decltype(*first + joiner) ret;

        if (first == last)
            return ret;
        
        ret += *first;
        for(++first; first != last; ++first) {
            ret += joiner;
            ret += *first;
        };
        return ret;
    }

    ARGUM_MOD_EXPORTED
    template<class It1, class It2, class Joiner, class Projection>
    auto join(It1 first, It2 last, Joiner joiner, Projection proj) -> decltype(proj(*first) + joiner) {
        
        decltype(proj(*first) + joiner) ret;

        if (first == last)
            return ret;
        
        ret += proj(*first);
        for(++first; first != last; ++first) {
            ret += joiner;
            ret += proj(*first);
        };
        return ret;
    }

    template<class Char>
    auto matchPrefix(std::basic_string_view<Char> value, std::basic_string_view<Char> prefix) -> bool {

        if (value.size() < prefix.size())
            return false;

        auto prefixCurrent = prefix.begin();
        auto prefixLast = prefix.end();
        auto valueCurrent = value.begin();
        for ( ; ; ) {
            if (*prefixCurrent != *valueCurrent)
                return false;
            if (++prefixCurrent == prefixLast)
                return true;
            ++valueCurrent;
        }
    }

    template<class Char>
    auto matchStrictPrefix(std::basic_string_view<Char> value, std::basic_string_view<Char> prefix) -> bool {

        if (value.size() <= prefix.size())
            return false;

        auto prefixCurrent = prefix.begin();
        auto prefixLast = prefix.end();
        auto valueCurrent = value.begin();
        for ( ; ; ) {
            if (*prefixCurrent != *valueCurrent)
                return false;
            if (++prefixCurrent == prefixLast)
                return true;
            ++valueCurrent;
        }
    }
    
}

#endif


namespace Argum {

    template<Character Char> struct CharConstants;

    #define ARGUM_DEFINE_CHAR_CONSTANTS(type, prefix) template<> struct CharConstants<type> { \
        static constexpr auto dash                          = prefix ## '-'; \
        static constexpr auto doubleDash                    = prefix ## "--"; \
        static constexpr auto assignment                    = prefix ## '='; \
        static constexpr auto slash                         = prefix ## '/'; \
        static constexpr auto colon                         = prefix ## ':'; \
        static constexpr auto braceOpen                     = prefix ## '{'; \
        static constexpr auto braceClose                    = prefix ## '}'; \
        static constexpr auto space                         = prefix ## ' '; \
        static constexpr auto squareBracketOpen             = prefix ## '['; \
        static constexpr auto squareBracketClose            = prefix ## ']'; \
        static constexpr auto pipe                          = prefix ## '|'; \
        static constexpr auto endl                          = prefix ## '\n'; \
        static constexpr auto ellipsis                      = prefix ## "..."; \
        static constexpr auto esc                           = prefix ## '\x1B'; \
        static constexpr auto semicolon                     = prefix ## ';'; \
        static constexpr auto letter_m                      = prefix ## 'm'; \
        static constexpr auto digit_0                       = prefix ## '0'; \
        static constexpr auto mustEscapeInRegex             = prefix ## R"([.^$|()\[\]{}*+?\\])"; \
        static constexpr auto regexEscapeReplacement        = prefix ## R"(\\&)"; \
        static constexpr auto falseNames                    = {prefix ## "0", prefix ## "false", prefix ## "off", prefix ## "no"}; \
        static constexpr auto trueNames                     = {prefix ## "1", prefix ## "true", prefix ## "on", prefix ## "yes"}; \
        \
        static auto isSpace(type c) -> bool; \
        static auto toLong(const type * str, type ** str_end, int base) -> long; \
        static auto toLongLong(const type * str, type ** str_end, int base) -> long long; \
        static auto toULong(const type * str, type ** str_end, int base) -> unsigned long; \
        static auto toULongLong(const type * str, type ** str_end, int base) -> unsigned long long; \
        static auto toFloat(const type * str, type ** str_end) -> float; \
        static auto toDouble(const type * str, type ** str_end) -> double; \
        static auto toLongDouble(const type * str, type ** str_end) -> long double; \
    };

    ARGUM_DEFINE_CHAR_CONSTANTS(char, )
    ARGUM_DEFINE_CHAR_CONSTANTS(wchar_t, L)

    #undef ARGUM_DEFINE_CHAR_CONSTANTS
    
    inline auto CharConstants<char>::isSpace(char c) -> bool { return isspace(c); }
    inline auto CharConstants<wchar_t>::isSpace(wchar_t c) -> bool { return iswspace(c); }

    inline auto CharConstants<char>::toLong(const char * str, char ** str_end, int base) -> long {
        return strtol(str, str_end, base);
    }
    inline auto CharConstants<wchar_t>::toLong(const wchar_t * str, wchar_t ** str_end, int base) -> long {
        return wcstol(str, str_end, base);
    }
    
    inline auto CharConstants<char>::toLongLong(const char * str, char ** str_end, int base) -> long long {
        return strtoll(str, str_end, base);
    }
    inline auto CharConstants<wchar_t>::toLongLong(const wchar_t * str, wchar_t ** str_end, int base) -> long long {
        return wcstoll(str, str_end, base);
    }

    inline auto CharConstants<char>::toULong(const char * str, char ** str_end, int base) -> unsigned long {
        return strtoul(str, str_end, base);
    }
    inline auto CharConstants<wchar_t>::toULong(const wchar_t * str, wchar_t ** str_end, int base) -> unsigned long {
        return wcstoul(str, str_end, base);
    }
    
    inline auto CharConstants<char>::toULongLong(const char * str, char ** str_end, int base) -> unsigned long long {
        return strtoull(str, str_end, base);
    }
    inline auto CharConstants<wchar_t>::toULongLong(const wchar_t * str, wchar_t ** str_end, int base) -> unsigned long long {
        return wcstoull(str, str_end, base);
    }

    inline auto CharConstants<char>::toFloat(const char * str, char ** str_end) -> float {
        return strtof(str, str_end);
    }
    inline auto CharConstants<wchar_t>::toFloat(const wchar_t * str, wchar_t ** str_end) -> float {
        return wcstof(str, str_end);
    }

    inline auto CharConstants<char>::toDouble(const char * str, char ** str_end) -> double {
        return strtod(str, str_end);
    }
    inline auto CharConstants<wchar_t>::toDouble(const wchar_t * str, wchar_t ** str_end) -> double {
        return wcstod(str, str_end);
    }
    
    inline auto CharConstants<char>::toLongDouble(const char * str, char ** str_end) -> long double {
        return strtold(str, str_end);
    }
    inline auto CharConstants<wchar_t>::toLongDouble(const wchar_t * str, wchar_t ** str_end) -> long double {
        return wcstold(str, str_end);
    }

    template<class Char>
    auto trimInPlace(std::basic_string<Char> & str) -> std::basic_string<Char> & {
        auto firstNotSpace = std::find_if(str.begin(), str.end(), [](const auto c) {
            return !CharConstants<Char>::isSpace(c);
        });
        str.erase(str.begin(), firstNotSpace);
        auto lastNotSpace = std::find_if(str.rbegin(), str.rend(), [](const auto c) {
            return !CharConstants<Char>::isSpace(c);
        }).base();
        str.erase(lastNotSpace, str.end());
        return str;
    }

}

#endif
#ifndef HEADER_ARGUM_MESSAGES_H_INCLUDED
#define HEADER_ARGUM_MESSAGES_H_INCLUDED




namespace Argum {

    template<Character Char> struct Messages;

    #if !ARGUM_OVERRIDE_MESSAGES

    #define ARGUM_DEFINE_MESSAGES(type, pr) template<> struct Messages<type> { \
        static constexpr auto unrecognizedOptionError()     { return pr ## "unrecognized option: {1}"; }\
        static constexpr auto ambiguousOptionError()        { return pr ## "ambigous option: {1}, candidates: {2}"; }\
        static constexpr auto missingOptionArgumentError()  { return pr ## "argument required for option: {1}"; }\
        static constexpr auto extraOptionArgumentError()    { return pr ## "extraneous argument for option: {1}"; }\
        static constexpr auto extraPositionalError()        { return pr ## "unexpected argument: {1}"; }\
        static constexpr auto validationError()             { return pr ## "invalid arguments: {1}"; }\
        static constexpr auto errorReadingResponseFile()    { return pr ## "error reading response file \"{1}\": {2}"; }\
        static constexpr auto option()                      { return pr ## "option"; }\
        static constexpr auto positionalArg()               { return pr ## "positional argument"; }\
        static constexpr auto usageStart()                  { return pr ## "Usage: "; }\
        static constexpr auto defaultArgName()              { return pr ## "ARG"; }\
        static constexpr auto positionalHeader()            { return pr ## "positional arguments:"; }\
        static constexpr auto optionsHeader()               { return pr ## "options:"; }\
        static constexpr auto listJoiner()                  { return pr ## ", "; }\
        static constexpr auto expectedValueNotPresent()     { return pr ## "expected value is not present"; }\
        static constexpr auto notANumber()                  { return pr ## "value \"{1}\" is not a number"; } \
        static constexpr auto outOfRange()                  { return pr ## "value \"{1}\" is out of range"; } \
        static constexpr auto notAValidChoice()             { return pr ## "value \"{1}\" is not one of the valid choices {{{2}}"; } \
        static constexpr auto itemUnrestricted()            { return pr ## "{1} {2} can occur any number of times"; }\
        static constexpr auto itemMustBePresent()           { return pr ## "{1} {2} must be present"; }\
        static constexpr auto itemMustNotBePresent()        { return pr ## "{1} {2} must not be present"; }\
        static constexpr auto itemMustBePresentGE(unsigned val) { \
            switch(val) { \
                case 0:   return itemUnrestricted(); \
                case 1:   return itemMustBePresent(); \
                default:  return pr ## "{1} {2} must occur at least {3} times"; \
            } \
        }\
        static constexpr auto itemMustBePresentG(unsigned val) { \
            switch(val) { \
                case 0:   return itemMustBePresent(); \
                default:  return pr ## "{1} {2} must occur more than {3} times"; \
            }\
        }\
        static constexpr auto itemMustBePresentLE(unsigned val) { \
            constexpr auto inf = std::numeric_limits<unsigned>::max();\
            switch(val) { \
                case 0:   return itemMustNotBePresent(); \
                case 1:   return pr ## "{1} {2} must occur at most once"; \
                default:  return pr ## "{1} {2} must occur at most {3} times"; \
                case inf: return itemUnrestricted(); \
            }\
        }\
        static constexpr auto itemMustBePresentL(unsigned val) { \
            constexpr auto inf = std::numeric_limits<unsigned>::max();\
            switch(val) { \
                case 0:   return itemMustNotBePresent();\
                case 1:   return itemMustNotBePresent();\
                default:  return pr ## "{1} {2} must occur less than {3} times";\
                case inf: return itemUnrestricted();\
            }\
        }\
        static constexpr auto itemMustBePresentEQ(unsigned val) { \
            switch(val) { \
                case 0:   return itemMustNotBePresent();\
                case 1:   return itemMustBePresent();\
                default:  return pr ## "{1} {2} must occur {3} times";\
            }\
        }\
        static constexpr auto itemMustBePresentNEQ(unsigned val) { \
            constexpr auto inf = std::numeric_limits<unsigned>::max();\
            switch(val) { \
                case 0:   return itemMustBePresent();\
                default: return pr ## "{1} {2} must NOT occur {3} times"; \
                case inf: return itemUnrestricted();\
            }\
        }\
    };

    ARGUM_DEFINE_MESSAGES(char, )
    ARGUM_DEFINE_MESSAGES(wchar_t, L)

    #undef ARGUM_DEFINE_MESSAGES

    #endif

}


#endif
#ifndef HEADER_ARGUM_FORMATTING_H_INCLUDED
#define HEADER_ARGUM_FORMATTING_H_INCLUDED


#ifndef HEADER_ARGUM_WCWIDTH_H_INCLUDED
#define HEADER_ARGUM_WCWIDTH_H_INCLUDED


/*
* This is an implementation of wcwidth() and wcswidth() (defined in
* IEEE Std 1002.1-2001) for Unicode.
*
* http://www.opengroup.org/onlinepubs/007904975/functions/wcwidth.html
* http://www.opengroup.org/onlinepubs/007904975/functions/wcswidth.html
*
* In fixed-width output devices, Latin characters all occupy a single
* "cell" position of equal width, whereas ideographic CJK characters
* occupy two such cells. Interoperability between terminal-line
* applications and (teletype-style) character terminals using the
* UTF-8 encoding requires agreement on which character should advance
* the cursor by how many cell positions. No established formal
* standards exist at present on which Unicode character shall occupy
* how many cell positions on character terminals. These routines are
* a first attempt of defining such behavior based on simple rules
* applied to data provided by the Unicode Consortium.
*
* For some graphical characters, the Unicode standard explicitly
* defines a character-cell width via the definition of the East Asian
* FullWidth (F), Wide (W), Half-width (H), and Narrow (Na) classes.
* In all these cases, there is no ambiguity about which width a
* terminal shall use. For characters in the East Asian Ambiguous (A)
* class, the width choice depends purely on a preference of backward
* compatibility with either historic CJK or Western practice.
* Choosing single-width for these characters is easy to justify as
* the appropriate long-term solution, as the CJK practice of
* displaying these characters as double-width comes from historic
* implementation simplicity (8-bit encoded characters were displayed
* single-width and 16-bit ones double-width, even for Greek,
* Cyrillic, etc.) and not any typographic considerations.
*
* Much less clear is the choice of width for the Not East Asian
* (Neutral) class. Existing practice does not dictate a width for any
* of these characters. It would nevertheless make sense
* typographically to allocate two character cells to characters such
* as for instance EM SPACE or VOLUME INTEGRAL, which cannot be
* represented adequately with a single-width glyph. The following
* routines at present merely assign a single-cell width to all
* neutral characters, in the interest of simplicity. This is not
* entirely satisfactory and should be reconsidered before
* establishing a formal standard in this area. At the moment, the
* decision which Not East Asian (Neutral) characters should be
* represented by double-width glyphs cannot yet be answered by
* applying a simple rule from the Unicode database content. Setting
* up a proper standard for the behavior of UTF-8 character terminals
* will require a careful analysis not only of each Unicode character,
* but also of each presentation form, something the author of these
* routines has avoided to do so far.
*
* http://www.unicode.org/unicode/reports/tr11/
*
* Markus Kuhn -- 2007-05-26 (Unicode 5.0)
*
* Permission to use, copy, modify, and distribute this software
* for any purpose and without fee is hereby granted. The author
* disclaims all warranties with regard to this software.
*
* Latest version: http://www.cl.cam.ac.uk/~mgk25/ucs/wcwidth.c
*/


namespace Argum::Wcwidth {

struct interval {
    int first;
    int last;
};

/* auxiliary function for binary search in interval table */
inline int bisearch(char32_t ucs, const struct interval *table, int max) {
    int min = 0;
    int mid;

    if (int(ucs) < table[0].first || int(ucs) > table[max].last)
        return 0;
    while (max >= min) {
        mid = (min + max) / 2;
        if (int(ucs) > table[mid].last)
            min = mid + 1;
        else if (int(ucs) < table[mid].first)
            max = mid - 1;
        else
            return 1;
    }

    return 0;
}


/* The following two functions define the column width of an ISO 10646
* character as follows:
*
*    - The null character (U+0000) has a column width of 0.
*
*    - Other C0/C1 control characters and DEL will lead to a return
*      value of -1.
*
*    - Non-spacing and enclosing combining characters (general
*      category code Mn or Me in the Unicode database) have a
*      column width of 0.
*
*    - SOFT HYPHEN (U+00AD) has a column width of 1.
*
*    - Other format characters (general category code Cf in the Unicode
*      database) and ZERO WIDTH SPACE (U+200B) have a column width of 0.
*
*    - Hangul Jamo medial vowels and final consonants (U+1160-U+11FF)
*      have a column width of 0.
*
*    - Spacing characters in the East Asian Wide (W) or East Asian
*      Full-width (F) category as defined in Unicode Technical
*      Report #11 have a column width of 2.
*
*    - All remaining characters (including all printable
*      ISO 8859-1 and WGL4 characters, Unicode control characters,
*      etc.) have a column width of 1.
*
* This implementation assumes that wchar_t characters are encoded
* in ISO 10646.
*/

inline int mk_wcwidth(char32_t ucs)
{
    /* sorted list of non-overlapping intervals of non-spacing characters */
    /* generated by "uniset +cat=Me +cat=Mn +cat=Cf -00AD +1160-11FF +200B c" */
    static constexpr struct interval combining[] = {
        { 0x0300, 0x036F }, { 0x0483, 0x0486 }, { 0x0488, 0x0489 },
        { 0x0591, 0x05BD }, { 0x05BF, 0x05BF }, { 0x05C1, 0x05C2 },
        { 0x05C4, 0x05C5 }, { 0x05C7, 0x05C7 }, { 0x0600, 0x0603 },
        { 0x0610, 0x0615 }, { 0x064B, 0x065E }, { 0x0670, 0x0670 },
        { 0x06D6, 0x06E4 }, { 0x06E7, 0x06E8 }, { 0x06EA, 0x06ED },
        { 0x070F, 0x070F }, { 0x0711, 0x0711 }, { 0x0730, 0x074A },
        { 0x07A6, 0x07B0 }, { 0x07EB, 0x07F3 }, { 0x0901, 0x0902 },
        { 0x093C, 0x093C }, { 0x0941, 0x0948 }, { 0x094D, 0x094D },
        { 0x0951, 0x0954 }, { 0x0962, 0x0963 }, { 0x0981, 0x0981 },
        { 0x09BC, 0x09BC }, { 0x09C1, 0x09C4 }, { 0x09CD, 0x09CD },
        { 0x09E2, 0x09E3 }, { 0x0A01, 0x0A02 }, { 0x0A3C, 0x0A3C },
        { 0x0A41, 0x0A42 }, { 0x0A47, 0x0A48 }, { 0x0A4B, 0x0A4D },
        { 0x0A70, 0x0A71 }, { 0x0A81, 0x0A82 }, { 0x0ABC, 0x0ABC },
        { 0x0AC1, 0x0AC5 }, { 0x0AC7, 0x0AC8 }, { 0x0ACD, 0x0ACD },
        { 0x0AE2, 0x0AE3 }, { 0x0B01, 0x0B01 }, { 0x0B3C, 0x0B3C },
        { 0x0B3F, 0x0B3F }, { 0x0B41, 0x0B43 }, { 0x0B4D, 0x0B4D },
        { 0x0B56, 0x0B56 }, { 0x0B82, 0x0B82 }, { 0x0BC0, 0x0BC0 },
        { 0x0BCD, 0x0BCD }, { 0x0C3E, 0x0C40 }, { 0x0C46, 0x0C48 },
        { 0x0C4A, 0x0C4D }, { 0x0C55, 0x0C56 }, { 0x0CBC, 0x0CBC },
        { 0x0CBF, 0x0CBF }, { 0x0CC6, 0x0CC6 }, { 0x0CCC, 0x0CCD },
        { 0x0CE2, 0x0CE3 }, { 0x0D41, 0x0D43 }, { 0x0D4D, 0x0D4D },
        { 0x0DCA, 0x0DCA }, { 0x0DD2, 0x0DD4 }, { 0x0DD6, 0x0DD6 },
        { 0x0E31, 0x0E31 }, { 0x0E34, 0x0E3A }, { 0x0E47, 0x0E4E },
        { 0x0EB1, 0x0EB1 }, { 0x0EB4, 0x0EB9 }, { 0x0EBB, 0x0EBC },
        { 0x0EC8, 0x0ECD }, { 0x0F18, 0x0F19 }, { 0x0F35, 0x0F35 },
        { 0x0F37, 0x0F37 }, { 0x0F39, 0x0F39 }, { 0x0F71, 0x0F7E },
        { 0x0F80, 0x0F84 }, { 0x0F86, 0x0F87 }, { 0x0F90, 0x0F97 },
        { 0x0F99, 0x0FBC }, { 0x0FC6, 0x0FC6 }, { 0x102D, 0x1030 },
        { 0x1032, 0x1032 }, { 0x1036, 0x1037 }, { 0x1039, 0x1039 },
        { 0x1058, 0x1059 }, { 0x1160, 0x11FF }, { 0x135F, 0x135F },
        { 0x1712, 0x1714 }, { 0x1732, 0x1734 }, { 0x1752, 0x1753 },
        { 0x1772, 0x1773 }, { 0x17B4, 0x17B5 }, { 0x17B7, 0x17BD },
        { 0x17C6, 0x17C6 }, { 0x17C9, 0x17D3 }, { 0x17DD, 0x17DD },
        { 0x180B, 0x180D }, { 0x18A9, 0x18A9 }, { 0x1920, 0x1922 },
        { 0x1927, 0x1928 }, { 0x1932, 0x1932 }, { 0x1939, 0x193B },
        { 0x1A17, 0x1A18 }, { 0x1B00, 0x1B03 }, { 0x1B34, 0x1B34 },
        { 0x1B36, 0x1B3A }, { 0x1B3C, 0x1B3C }, { 0x1B42, 0x1B42 },
        { 0x1B6B, 0x1B73 }, { 0x1DC0, 0x1DCA }, { 0x1DFE, 0x1DFF },
        { 0x200B, 0x200F }, { 0x202A, 0x202E }, { 0x2060, 0x2063 },
        { 0x206A, 0x206F }, { 0x20D0, 0x20EF }, { 0x302A, 0x302F },
        { 0x3099, 0x309A }, { 0xA806, 0xA806 }, { 0xA80B, 0xA80B },
        { 0xA825, 0xA826 }, { 0xFB1E, 0xFB1E }, { 0xFE00, 0xFE0F },
        { 0xFE20, 0xFE23 }, { 0xFEFF, 0xFEFF }, { 0xFFF9, 0xFFFB },
        { 0x10A01, 0x10A03 }, { 0x10A05, 0x10A06 }, { 0x10A0C, 0x10A0F },
        { 0x10A38, 0x10A3A }, { 0x10A3F, 0x10A3F }, { 0x1D167, 0x1D169 },
        { 0x1D173, 0x1D182 }, { 0x1D185, 0x1D18B }, { 0x1D1AA, 0x1D1AD },
        { 0x1D242, 0x1D244 }, { 0xE0001, 0xE0001 }, { 0xE0020, 0xE007F },
        { 0xE0100, 0xE01EF }
    };

    /* test for 8-bit control characters */
    if (ucs == 0)
        return 0;
    if (ucs < 32 || (ucs >= 0x7f && ucs < 0xa0))
        return -1;

    /* binary search in table of non-spacing characters */
    if (bisearch(ucs, combining,
        sizeof(combining) / sizeof(struct interval) - 1))
        return 0;

    /* if we arrive here, ucs is not a combining or C0/C1 control character */

    return 1 + 
        (ucs >= 0x1100 &&
            (ucs <= 0x115f ||                    /* Hangul Jamo init. consonants */
                ucs == 0x2329 || ucs == 0x232a ||
                (ucs >= 0x2e80 && ucs <= 0xa4cf &&
                    ucs != 0x303f) ||                  /* CJK ... Yi */
                (ucs >= 0xac00 && ucs <= 0xd7a3) || /* Hangul Syllables */
                (ucs >= 0xf900 && ucs <= 0xfaff) || /* CJK Compatibility Ideographs */
                (ucs >= 0xfe10 && ucs <= 0xfe19) || /* Vertical forms */
                (ucs >= 0xfe30 && ucs <= 0xfe6f) || /* CJK Compatibility Forms */
                (ucs >= 0xff00 && ucs <= 0xff60) || /* Fullwidth Forms */
                (ucs >= 0xffe0 && ucs <= 0xffe6) ||
                (ucs >= 0x20000 && ucs <= 0x2fffd) ||
                (ucs >= 0x30000 && ucs <= 0x3fffd)));
}


inline int mk_wcswidth(const char32_t *pwcs, size_t n)
{
    int w, width = 0;

    for (;*pwcs && n-- > 0; pwcs++)
        if ((w = mk_wcwidth(*pwcs)) < 0)
            return -1;
        else
            width += w;

    return width;
}

inline int wcswidth(const wchar_t *pwcs, size_t n) {

    if constexpr (sizeof(wchar_t) == 2) {

        int width = 0;
        uint32_t acc = 0;
        bool need_trail = false;
        for (;*pwcs && n-- > 0; pwcs++) {
            auto c = uint16_t(*pwcs);

            if (need_trail) {
                if ((c & 0xFC00) == 0xDC00) {
                    need_trail = false;
                    acc = ((acc - 0xD800) << 10) + (uint32_t(c) - 0xDC00) + 0x0010000;
                } else {
                    return -1;
                }
            } else if ((c & 0xFC00) == 0xD800) {
                acc = c;
                need_trail = true;
                continue;
            } else {
                acc = c;
            }

            if (int w = mk_wcwidth(char32_t(acc)); w < 0)
                return -1;
            else
                width += w;
        }

        if (need_trail)
            return -1;

        return width;
        
    } else if constexpr (sizeof(wchar_t) == 4) {
        
        return mk_wcswidth((const char32_t *)pwcs, n);
    } else {
        return -1;
    }

}

}

#endif 



namespace Argum {

    template<class T>
    auto append(std::string & dest, T val) -> std::string &
    requires(std::is_same_v<decltype(std::to_string(val)), std::string>) {
        return dest += std::to_string(val);
    }

    template<class T>
    auto append(std::wstring & dest, T val) -> std::wstring &
    requires(std::is_same_v<decltype(std::to_wstring(val)), std::wstring>)  {
        return dest += std::to_wstring(val);
    }

    template<Character Char, class T>
    auto append(std::basic_string<Char> & dest, T && val) -> std::basic_string<Char> &
    requires(std::is_same_v<decltype(dest += std::forward<T>(val)), decltype(dest)>) {
        return dest += std::forward<T>(val);
    }

    inline auto append(std::string & dest, const wchar_t * wstr) -> std::string & {
        
        mbstate_t state = mbstate_t();
        size_t len = wcsrtombs(nullptr, &wstr, 0, &state);
        if (len == size_t(-1))
            return dest;
        auto currentSize = dest.size();
        dest.resize(currentSize + len + 1);
        wcsrtombs(&dest[currentSize], &wstr, len + 1, &state);
        dest.resize(currentSize + len);
        return dest;
    }

    inline auto append(std::string & dest, const std::wstring & str) -> std::string & {
        return append(dest, str.c_str());
    }

    inline auto append(std::string & dest, const std::wstring_view & str) -> std::string & {
        return append(dest, std::wstring(str));
    }

    inline auto append(std::wstring & dest, const char * str) -> std::wstring & {
        
        mbstate_t state = mbstate_t();
        size_t len = mbsrtowcs(nullptr, &str, 0, &state);
        if (len == size_t(-1))
            return dest;
        auto currentSize = dest.size();
        dest.resize(currentSize + len + 1);
        mbsrtowcs(&dest[currentSize], &str, len + 1, &state);
        dest.resize(currentSize + len);
        return dest;
    }

    inline auto append(std::wstring & dest, const std::string & str) -> std::wstring & {
        return append(dest, str.c_str());
    }

    inline auto append(std::wstring & dest, const std::string_view & str) -> std::wstring & {
        return append(dest, std::string(str));
    }

    template<class T, class Char>
    concept StringAppendable = Character<Char> && 
        requires(std::basic_string<Char> & str, T && val) {
            { append(str, std::forward<T>(val)) } -> std::same_as<std::basic_string<Char> &>;
    };

    namespace Impl {
        template <typename T>
        concept has_wcwidth = requires(T * t) {
            { wcswidth(t, size_t(0)) };
        };

        constexpr bool wcswidthPresent = has_wcwidth<wchar_t>;

        template<class T = wchar_t>
        inline auto simpleWidth(const T * str, size_t size) -> int {
            if constexpr (wcswidthPresent) {
                return wcswidth(str, size);
            } else {
                return Argum::Wcwidth::wcswidth(str, size);
            }
        }
    }

    inline auto stringWidth(const std::wstring_view & str) -> unsigned {

        int res = Impl::simpleWidth(str.data(), str.size());
        if (res >= 0)
            return res;
        
        std::wstring stripped;
        enum {
            stateNormal,
            stateEsc,
            stateControlStart,
            stateControlIntermediate
        } state = stateNormal;
        for (wchar_t c: str) {
            switch(state) {
            break; case stateNormal: restart:
                if (c == L'\x1b') {
                    state = stateEsc;
                    continue;
                }
                stripped += c;

            break; case stateEsc:
                if (c == L'[') {
                    state = stateControlStart;
                    continue;
                }
                state = stateNormal;
                goto restart;

            break; case stateControlStart:
                if (c >= 0x30 && c <= 0x3F) {
                    continue;
                }
                if (c >= 0x20 && c <= 0x2F) {
                    state = stateControlIntermediate;
                    continue;
                }
                if (c >= 0x40 && c <= 0x7E) {
                    state = stateNormal;
                    continue;
                }
                state = stateNormal;
                goto restart;

            break; case stateControlIntermediate:
                if (c >= 0x20 && c <= 0x2F) {
                    state = stateControlIntermediate;
                    continue;
                }
                if (c >= 0x40 && c <= 0x7E) {
                    state = stateNormal;
                    continue;
                }
                state = stateNormal;
                goto restart;
            }
        }

        res = Impl::simpleWidth(stripped.data(), stripped.size());
        if (res >= 0)
            return res;
        
        return unsigned(stripped.size());
    }

    inline auto stringWidth(const std::string_view & str) -> unsigned {

        std::wstring sigh;
        append(sigh, str);
        return stringWidth(sigh);
    }

    template<class Char, StringAppendable<Char> T>
    auto toString(T && val) -> std::basic_string<Char> {
        std::basic_string<Char> ret;
        append(ret, std::forward<T>(val));
        return ret;
    }
    
    template<Character Char>
    auto appendByIndex(std::basic_string<Char> & /*dest*/, size_t /*idx*/) {

        //do nothing
    }

    template<Character Char, StringAppendable<Char> First, StringAppendable<Char>... Rest>
    auto appendByIndex(std::basic_string<Char> & dest, size_t idx, First && first, Rest && ...rest) {

        if (idx == 0) {
            append(dest, std::forward<First>(first));
            return;
        }

        appendByIndex(dest, idx - 1, rest...);
    }

    template<Character Char>
    auto parseFormatPlaceholder(const Char * first, const Char * last) -> std::optional<size_t> {
            
        [[maybe_unused]] char buffer[16 * MB_LEN_MAX]; //16 characters for placeholder number ought to be enough :)
        const char * convertFirst;
        const char * convertLast;
        
        if constexpr (!std::is_same_v<Char, char>) {

            char * const bufferStart = std::begin(buffer);
            char * const bufferEnd = std::end(buffer);
            
            char * bufferCurrent = bufferStart;
            std::mbstate_t state = {};
            for ( ; first != last; ++first) {
                if (size_t(bufferEnd - bufferCurrent) < size_t(MB_CUR_MAX))
                    return std::nullopt;

                size_t written;
                if constexpr (std::is_same_v<Char, wchar_t>)
                    written = wcrtomb(bufferCurrent, *first, &state);
                else if constexpr (std::is_same_v<Char, char8_t>)
                    written = c8rtomb(bufferCurrent, *first, &state);
                else if constexpr (std::is_same_v<Char, char16_t>)
                    written = c16rtomb(bufferCurrent, *first, &state);
                else if constexpr (std::is_same_v<Char, char32_t>)
                    written = c32rtomb(bufferCurrent, *first, &state);

                if (written == size_t(-1))
                    return std::nullopt;
                bufferCurrent += written;
            }
            convertFirst = bufferStart;
            convertLast = bufferCurrent;

        } else {
            
            convertFirst = first;
            convertLast = last;
        }
        
        size_t ret;
        auto res = std::from_chars(convertFirst, convertLast, ret);
        if (res.ptr != convertLast || res.ec != std::errc())
            return std::nullopt;
        return ret;
    }

    template<StringLike Fmt, StringAppendable<CharTypeOf<Fmt>>... Args>
    auto format(Fmt && fmt, Args && ...args)  {

        using CharType = CharTypeOf<Fmt>;
        using StringViewType = std::basic_string_view<CharType>;
        using StringType = std::basic_string<CharType>;
        using CharConstants = CharConstants<CharType>;

        StringViewType fmtView(fmt);

        StringType ret;

        auto outStart = fmtView.data();
        auto current = outStart;
        const auto last = current + fmtView.size();

        auto flush = [&](auto end) {
            assert(end >= outStart);
            append(ret, StringViewType(outStart, size_t(end - outStart)));
        };

        while(current != last) {

            if (*current == CharConstants::braceOpen) {
                    ++current; 
                    if (current == last) 
                        break;
                    if (*current == CharConstants::braceOpen) {
                        flush(current);
                        ++current;
                        outStart = current;
                        continue;
                    }

                    auto placeholderEnd = std::find(current, last, CharConstants::braceClose);
                    if (placeholderEnd == last) 
                        continue;
                    auto maybeArgIdx = parseFormatPlaceholder(current, placeholderEnd);
                    if (!maybeArgIdx || *maybeArgIdx == 0 || *maybeArgIdx > sizeof...(Args)) 
                        continue;

                    flush(current - 1);
                    appendByIndex(ret, *maybeArgIdx - 1, args...);
                    current = placeholderEnd + 1;
                    outStart = current;
            } else {
                ++current; 
            }
        }
        flush(last);
        return ret;
    }

    template<StringLike T>
    auto wordWrap(T && input, unsigned maxLength, 
                  unsigned indent = 0, unsigned firstLineOffset = 0) -> std::basic_string<CharTypeOf<T>> {

        using Char = CharTypeOf<T>;

        std::basic_string_view<Char> str(std::forward<T>(input));
        
        if (maxLength == 0 || str.empty())
            return {};

        if (indent >= maxLength)
            indent = maxLength - 1;

        std::basic_string<Char> ret;
        ret.reserve(str.size());
        std::basic_string<Char> prefix;
        bool firstLine = true;

        unsigned curMaxLen;
        if (firstLineOffset >= maxLength) {
            ret += CharConstants<Char>::endl;
            prefix = std::basic_string<Char>(indent, CharConstants<char>::space);
            ret += prefix;
            curMaxLen = maxLength - indent;
            firstLine = false;
        } else {
            curMaxLen = maxLength - firstLineOffset;
        }
        
        for ( ; ; ) {
            auto eolPos = str.find(CharConstants<Char>::endl);
            auto line = str.substr(0, eolPos);
            bool needLineBreak = (eolPos != str.npos);
            auto width = stringWidth(line);
            while (width > curMaxLen) {
                auto spacePos = line.rfind(CharConstants<Char>::space);
                if (spacePos == line.npos)
                    break;
                line = line.substr(0, spacePos);
                needLineBreak = true;
                width = stringWidth(line);
            }
            ret += prefix;
            ret += line;
            if (needLineBreak) {
                ret += CharConstants<Char>::endl;
                str = str.substr(line.size() + 1);
            } else {
                str = str.substr(line.size());
            }
            if (str.empty())
                break;

            if (firstLine) {
                prefix = std::basic_string<Char>(indent, CharConstants<char>::space);
                curMaxLen = maxLength - indent;
                firstLine = false;
            }
        }
        return ret;
    }

}

#endif
#ifndef HEADER_ARGUM_VALIDATORS_H_INCLUDED
#define HEADER_ARGUM_VALIDATORS_H_INCLUDED

#ifndef HEADER_ARGUM_FLAT_MAP_H_INCLUDED
#define HEADER_ARGUM_FLAT_MAP_H_INCLUDED



namespace Argum {

    template<class T>
    class FlatSet {
    public:
        using value_type = T;
        using iterator = typename std::vector<value_type>::iterator;
        using const_iterator = typename std::vector<value_type>::const_iterator;

    private:
        struct Comparator {
            template<class X>
            auto operator()(const value_type & lhs, const X & rhs) const {
                return lhs < rhs;
            }
            template<class X>
            auto operator()(const X & lhs, const value_type & rhs) const {
                return lhs < rhs;
            }
            auto operator()(const value_type & lhs, const value_type & rhs) const {
                return lhs < rhs;
            }
        };

    public:
        FlatSet() = default;
        
        FlatSet(const value_type & value): m_data(1, value) {
        }
        FlatSet(value_type && value) {
            this->m_data.emplace_back(std::move(value));
        }
        FlatSet(std::initializer_list<value_type> values) : m_data(values) {
            std::sort(this->m_data.begin(), this->m_data.end(), Comparator());
        }

        auto insert(T val) -> std::pair<iterator, bool> {
            auto it = std::lower_bound(this->m_data.begin(), this->m_data.end(), val, Comparator());
            if (it != this->m_data.end() && *it == val)
                return {std::move(it), false};
            it = this->m_data.emplace(it, std::move(val));
            return {std::move(it), true};
        }

        template<class Arg>
        auto find(const Arg & key) const -> const_iterator {
            auto it = std::lower_bound(this->m_data.begin(), this->m_data.end(), key, Comparator());
            if (it == this->m_data.end() || *it != key)
                return this->m_data.end();
            return it;
        }

        iterator begin() noexcept { return this->m_data.begin(); }
        const_iterator begin() const noexcept { return this->m_data.begin(); }
        const_iterator cbegin() const noexcept { return this->m_data.begin(); }
        iterator end() noexcept { return this->m_data.end(); }
        const_iterator end() const noexcept { return this->m_data.end(); }
        const_iterator cend() const noexcept { return this->m_data.end(); }

        size_t size() const { return m_data.size(); }
        size_t empty() const { return m_data.empty(); }

        void clear() { m_data.clear(); }

    private:
        std::vector<value_type> m_data;
    };

    template<class Key, class Value>
    class FlatMap {
    public:
        class value_type : private std::pair<Key, Value> {
            
        public:
            using std::pair<Key, Value>::pair;
            
            const Key & key() const noexcept {
                return this->first;
            }
            
            const Value & value() const noexcept {
                return this->second;
            }
            
            Value & value() noexcept {
                return this->second;
            }
        };
        using iterator = typename std::vector<value_type>::iterator;
        using const_iterator = typename std::vector<value_type>::const_iterator;
    private:
        struct Comparator {
            template<class X>
            auto operator()(const value_type & lhs, const X & rhs) const {
                return lhs.key() < rhs;
            }
            template<class X>
            auto operator()(const X & lhs, const value_type & rhs) const {
                return lhs < rhs.key();
            }
            auto operator()(const value_type & lhs, const value_type & rhs) {
                return lhs.key() < rhs.key();
            }
        };
    public:
        auto add(Key key, Value val) -> std::pair<iterator, bool> {
            auto it = std::lower_bound(this->m_data.begin(), this->m_data.end(), key, Comparator());
            if (it != this->m_data.end() && it->key() == key)
                return {std::move(it), false};
            it = this->m_data.emplace(it, std::move(key), std::move(val));
            return {std::move(it), true};
        }

        template<class KeyArg>
        auto operator[](const KeyArg & key) -> Value & {
            auto it = std::lower_bound(this->m_data.begin(), this->m_data.end(), key, Comparator());
            if (it == this->m_data.end() || it->key() != key)
                it = this->m_data.emplace(it, std::move(key), Value());
            return it->value();
        }

        auto valueByIndex(size_t idx) const -> Value & {
            return this->m_data[idx].value();
        }

        template<class KeyArg>
        auto find(const KeyArg & key) const -> const_iterator {
            auto it = std::lower_bound(this->m_data.begin(), this->m_data.end(), key, Comparator());
            if (it == this->m_data.end() || it->key() != key)
                return this->m_data.end();
            return it;
        }

        template<class KeyArg>
        auto lower_bound(const KeyArg & key) const -> const_iterator {
            return std::lower_bound(this->m_data.begin(), this->m_data.end(), key, Comparator());
        }

        iterator begin() noexcept { return this->m_data.begin(); }
        const_iterator begin() const noexcept { return this->m_data.begin(); }
        const_iterator cbegin() const noexcept { return this->m_data.begin(); }
        iterator end() noexcept { return this->m_data.end(); }
        const_iterator end() const noexcept { return this->m_data.end(); }
        const_iterator cend() const noexcept { return this->m_data.end(); }

        size_t size() const { return m_data.size(); }
        size_t empty() const { return m_data.empty(); }

        void clear() { m_data.clear(); }

    private:
        std::vector<value_type> m_data;
    };

    template<class Key, class Value, class Arg>
    auto findMatchOrMatchingPrefixRange(const FlatMap<Key, Value> & map, Arg arg) -> 
            std::pair<typename FlatMap<Key, Value>::const_iterator,
                      typename FlatMap<Key, Value>::const_iterator> {
        
        const auto end = map.end();

        auto first = map.lower_bound(arg);

        //at this point it points to element grater or equal than arg
        //so range [it, last) either:
        // 1. starts with an exact match
        // 2. contains a run of items which start with arg followed by ones which don't
        // 3. is empty

        //discard the empty range
        if (first == end)
            return {end, end};

        //if exact match, return it now
        if (first->key() == arg)
            return {first, first + 1};
    
        //figure out case 2 and the end of the prefix run
        auto last = first;
        do {
            auto compRes = last->key().compare(0, arg.size(), arg);
            assert(compRes >= 0); //invariant of range [it, last)
            
            //break if we reached the end of prefix run
            if (compRes != 0)
                break;
            
            ++last;
            
        } while (last != end);

        return {first, last};
    }

}


#endif 


namespace Argum {

    ARGUM_MOD_EXPORTED
    template<class Char>
    class BasicValidationData {
    public:
        auto optionCount(std::basic_string_view<Char> name) const -> unsigned {
            auto it = this->m_optionCounts.find(name);
            if (it == this->m_optionCounts.end())
                return 0;
            return it->value();
        }
        auto optionCount(std::basic_string_view<Char> name) -> unsigned & {
            return this->m_optionCounts[name];
        }
        auto positionalCount(std::basic_string_view<Char> name) const -> unsigned {
            auto it = this->m_positionalCounts.find(name);
            if (it == this->m_positionalCounts.end())
                return 0;
            return it->value();
        }
        auto positionalCount(std::basic_string_view<Char> name) -> unsigned & {
            return this->m_positionalCounts[name];
        }
    private:
        FlatMap<std::basic_string<Char>, unsigned> m_optionCounts;
        FlatMap<std::basic_string<Char>, unsigned> m_positionalCounts;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(ValidationData)

    template<class T, class Char>
    concept ParserValidator = std::is_invocable_r_v<bool, T, const BasicValidationData<Char>>;
    
    template<class T, class Char>
    concept DescribableParserValidator = ParserValidator<T, Char> && 
        requires(const T & val) {
            { describe(val) } -> std::same_as<std::basic_string<Char>>;
        };
    
    template<class Char, class First, class... Rest>
    constexpr bool CompatibleParserValidators = ParserValidator<First, Char> && (ParserValidator<Rest, Char> && ...);
    
    template<class T>
    concept AnyParserValidator =
        ParserValidator<T, char> ||
        ParserValidator<T, wchar_t>;
   
    
    template<class First, class... Rest>
    constexpr bool AnyCompatibleParserValidators =
        CompatibleParserValidators<char, First, Rest...> ||
        CompatibleParserValidators<wchar_t, First, Rest...>;

    //MARK: - NotValidator

    template<AnyParserValidator Impl>
    class NotValidator {
    public:
        template<class X>
        NotValidator(X && x) requires(std::is_convertible_v<X &&, Impl>): m_impl(std::forward<X>(x)) {
        }
        
        template<class Char>
        requires(ParserValidator<Impl, Char>)
        auto operator()(const BasicValidationData<Char> & data) const -> bool {
            return !this->m_impl(data);
        }
        
        auto operator!() const -> const Impl & {
            return m_impl;
        }
    private:
        Impl m_impl;
    };

    //MARK: - OppositeOf
    
    ARGUM_MOD_EXPORTED
    template<AnyParserValidator Validator>
    auto operator!(const Validator & val) {
        return NotValidator<std::remove_cvref_t<Validator>>(val);
    }

    ARGUM_MOD_EXPORTED
    template<AnyParserValidator Validator>
    auto oppositeOf(Validator && val) {
        return NotValidator<std::remove_cvref_t<Validator>>(val);
    }

    //MARK: - CombinedValidator

    enum class ValidatorCombination : int {
        And,
        Or,
        OnlyOne,
        OneOrNone,
        AllOrNone
    };

    template<ValidatorCombination Comb, class... Args>
    requires(sizeof...(Args) > 1 && AnyCompatibleParserValidators<Args...>)
    class CombinedValidator {
    public:
        using TupleType = std::tuple<Args...>;
    public:
        CombinedValidator(std::tuple<Args...> && args): m_items(std::move(args)) {
        }
        template<class... CArgs>
        requires(std::is_constructible_v<TupleType, CArgs && ...>)
        CombinedValidator(CArgs && ...args): m_items(std::forward<CArgs>(args)...) {  
        }

        template<class Char>
        requires(CompatibleParserValidators<Char, Args...>)
        auto operator()(const BasicValidationData<Char> & data) const -> bool {
            if constexpr (Comb == ValidatorCombination::And) {
                return std::apply([&data](const Args & ...args) -> bool { return (args(data) && ...); }, this->m_items);
            } else if constexpr (Comb == ValidatorCombination::Or) {
                return std::apply([&data](const Args & ...args) -> bool { return (args(data) || ...); }, this->m_items);
            } else if constexpr (Comb == ValidatorCombination::OnlyOne) {
                unsigned counter = 0;
                return this->evalOneOrNone<sizeof...(Args)>(data, counter) && counter == 1;
            } else if constexpr (Comb == ValidatorCombination::OneOrNone) {
                unsigned counter = 0;
                return this->evalOneOrNone<sizeof...(Args)>(data, counter);
            } else if constexpr (Comb == ValidatorCombination::AllOrNone) {
                unsigned counter = 0;
                return this->evalAllOrNone<sizeof...(Args)>(data, counter);
            }
        }

        auto operator!() const {
            
            if constexpr (Comb == ValidatorCombination::And || Comb == ValidatorCombination::Or)
                return std::apply([](const Args & ...args) {
                    if constexpr (Comb == ValidatorCombination::And)
                        return CombinedValidator<ValidatorCombination::Or, decltype(!args)...>(!args...);
                    else if constexpr (Comb == ValidatorCombination::Or)
                        return CombinedValidator<ValidatorCombination::And, decltype(!args)...>(!args...);
                }, this->m_items);
            else 
                return NotValidator<CombinedValidator>(*this);
        }

        auto items() const -> const TupleType & {
            return m_items;
        }

    private:

        template<size_t N, class Char>
        requires(CompatibleParserValidators<Char, Args...>)
        auto evalOneOrNone(const BasicValidationData<Char> & data, unsigned & counter) const {
            
            constexpr auto idx = sizeof...(Args) - N;

            if (std::get<idx>(this->m_items)(data)) {
                if (++counter > 1)
                    return false;
            }

            if constexpr (N > 1) {
                return this->evalOneOrNone<N - 1>(data, counter);
            } else {
                return true;
            }
        }

        template<size_t N, class Char>
        requires(CompatibleParserValidators<Char, Args...>)
        auto evalAllOrNone(const BasicValidationData<Char> & data, unsigned & counter) const {

            constexpr auto idx = sizeof...(Args) - N;
            
            counter += std::get<idx>(this->m_items)(data);
            
            if (counter != 0 && counter != idx + 1)
                return false;

            if constexpr (N > 1) {
                return this->evalAllOrNone<N - 1>(data, counter);
            } else {
                return true;
            }
        }

    private:
        TupleType m_items;
    };

    template<ValidatorCombination Comb, class V1, class V2>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<V1>, std::remove_cvref_t<V2>>)
    auto combine(V1 && v1, V2 && v2)  {
        using R1 = std::remove_cvref_t<V1>;
        using R2 = std::remove_cvref_t<V2>;
        return CombinedValidator<Comb, R1, R2>(std::forward<V1>(v1), std::forward<V2>(v2));
    }

    template<ValidatorCombination Comb, class V1, class... Args2>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<V1>, Args2...>)
    auto combine(V1 && v1, const CombinedValidator<Comb, Args2...> & v2)  {
        using R1 = std::remove_cvref_t<V1>;
        return CombinedValidator<Comb, R1, Args2...>(
            std::tuple_cat(std::tuple<R1>(std::forward<V1>(v1)), v2.items())
        );
    }
    
    template<ValidatorCombination Comb, class V1, class... Args2>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<V1>, Args2...>)
    auto combine(V1 && v1, CombinedValidator<Comb, Args2...> && v2)  {
        using R1 = std::remove_cvref_t<V1>;
        return CombinedValidator<Comb, R1, Args2...>(
            std::tuple_cat(std::tuple<R1>(std::forward<V1>(v1)), std::move(v2.items()))
        );
    }

    template<ValidatorCombination Comb, class V2, class... Args1>
    requires(AnyCompatibleParserValidators<Args1..., std::remove_cvref_t<V2>>)
    auto combine(const CombinedValidator<Comb, Args1...> & v1, V2 && v2)  {
        using R2 = std::remove_cvref_t<V2>;
        return CombinedValidator<Comb, Args1..., R2>(
            std::tuple_cat(v1.items(), std::tuple<R2>(std::forward<V2>(v2)))
        );
    }
    
    template<ValidatorCombination Comb, class V2, class... Args1>
    requires(AnyCompatibleParserValidators<Args1..., std::remove_cvref_t<V2>>)
    auto combine(CombinedValidator<Comb, Args1...> && v1, V2 && v2)  {
        using R2 = std::remove_cvref_t<V2>;
        return CombinedValidator<Comb, Args1..., R2>(
            std::tuple_cat(std::move(v1.items()), std::tuple<R2>(std::forward<V2>(v2)))
        );
    }

    template<ValidatorCombination Comb, class... Args1, class... Args2>
    requires(AnyCompatibleParserValidators<Args1..., Args2...>)
    auto combine(const CombinedValidator<Comb, Args1...> & v1, const CombinedValidator<Comb, Args2...> & v2)  {
        return CombinedValidator<Comb, Args1..., Args2...>(
            std::tuple_cat(v1.items(), v2.items())
        );
    }
    
    template<ValidatorCombination Comb, class... Args1, class... Args2>
    requires(AnyCompatibleParserValidators<Args1..., Args2...>)
    auto combine(CombinedValidator<Comb, Args1...> && v1, const CombinedValidator<Comb, Args2...> & v2)  {
        return CombinedValidator<Comb, Args1..., Args2...>(
            std::tuple_cat(std::move(v1.items()), v2.items())
        );
    }
    
    template<ValidatorCombination Comb, class... Args1, class... Args2>
    requires(AnyCompatibleParserValidators<Args1..., Args2...>)
    auto combine(const CombinedValidator<Comb, Args1...> & v1, CombinedValidator<Comb, Args2...> && v2)  {
        return CombinedValidator<Comb, Args1..., Args2...>(
            std::tuple_cat(v1.items(), std::move(v2.items()))
        );
    }
    
    template<ValidatorCombination Comb, class... Args1, class... Args2>
    requires(AnyCompatibleParserValidators<Args1..., Args2...>)
    auto combine(CombinedValidator<Comb, Args1...> && v1, CombinedValidator<Comb, Args2...> && v2)  {
        return CombinedValidator<Comb, Args1..., Args2...>(
            std::tuple_cat(std::move(v1.items()), std::move(v2.items()))
        );
    }

    
    //MARK: - AllOf

    ARGUM_MOD_EXPORTED
    template<class V1, class V2>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<V1>, std::remove_cvref_t<V2>>)
    auto operator&&(V1 && v1, V2 && v2)  {
        return combine<ValidatorCombination::And>(std::forward<V1>(v1), std::forward<V2>(v2));
    }

    ARGUM_MOD_EXPORTED
    template<class First, class... Rest>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>)
    auto allOf(First && first, Rest && ...rest)  {
        return (std::forward<First>(first) &&  ... && std::forward<Rest>(rest));
    }

    //MARK: - AnyOf

    ARGUM_MOD_EXPORTED
    template<class V1, class V2>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<V1>, std::remove_cvref_t<V2>>)
    auto operator||(V1 && v1, V2 && v2)  {
        return combine<ValidatorCombination::Or>(std::forward<V1>(v1), std::forward<V2>(v2));
    }

    ARGUM_MOD_EXPORTED
    template<class First, class... Rest>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>)
    auto anyOf(First && first, Rest && ...rest)  {
        return (std::forward<First>(first) || ... || std::forward<Rest>(rest));
    }

    //MARK: - NoneOf

    ARGUM_MOD_EXPORTED
    template<class First, class... Rest>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>)
    auto noneOf(First && first, Rest && ...rest)  {
        return allOf(!std::forward<First>(first), !std::forward<Rest>(rest)...);
    }

    //MARK: - OnlyOneOf

    ARGUM_MOD_EXPORTED
    template<class First, class... Rest>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>)
    auto onlyOneOf(First && first, Rest && ...rest)  {
        
        return CombinedValidator<ValidatorCombination::OnlyOne, std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>(std::forward<First>(first), std::forward<Rest>(rest)...);
    }

    //MARK: - OneOrNoneOf

    ARGUM_MOD_EXPORTED
    template<class First, class... Rest>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>)
    auto oneOrNoneOf(First && first, Rest && ...rest)  {
        return CombinedValidator<ValidatorCombination::OneOrNone, std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>(std::forward<First>(first), std::forward<Rest>(rest)...);
    }

    //MARK: - AllOrNoneOf

    ARGUM_MOD_EXPORTED
    template<class First, class... Rest>
    requires(AnyCompatibleParserValidators<std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>)
    auto allOrNoneOf(First && first, Rest && ...rest)  {
        return CombinedValidator<ValidatorCombination::AllOrNone, std::remove_cvref_t<First>, std::remove_cvref_t<Rest>...>(std::forward<First>(first), std::forward<Rest>(rest)...);
    }

    //MARK: - Occurence Validators

    template<class Char, bool IsOption, class Comp>
    requires(std::is_same_v<Comp, std::greater_equal<unsigned>> ||
             std::is_same_v<Comp, std::less_equal<unsigned>> || 
             std::is_same_v<Comp, std::greater<unsigned>> ||
             std::is_same_v<Comp, std::less<unsigned>> ||
             std::is_same_v<Comp, std::equal_to<unsigned>> ||
             std::is_same_v<Comp, std::not_equal_to<unsigned>>)
    class ItemOccurs {
    public:
        ItemOccurs(std::basic_string_view<Char> name, unsigned count) : m_name(name), m_count(count) {}

        auto operator()(const BasicValidationData<Char> & data) const -> bool {
            if constexpr (IsOption)
                return Comp()(data.optionCount(this->m_name), this->m_count);
            else
                return Comp()(data.positionalCount(this->m_name), this->m_count);
        }
        
        auto operator!() const {
            if constexpr (std::is_same_v<Comp, std::greater_equal<unsigned>>)
                return ItemOccurs<Char, IsOption, std::less<unsigned>>(this->m_name, this->m_count);
            else if constexpr (std::is_same_v<Comp, std::less_equal<unsigned>>)
                return ItemOccurs<Char, IsOption, std::greater<unsigned>>(this->m_name, this->m_count);
            else if constexpr (std::is_same_v<Comp, std::greater<unsigned>>)
                return ItemOccurs<Char, IsOption, std::less_equal<unsigned>>(this->m_name, this->m_count);
            else if constexpr (std::is_same_v<Comp, std::less<unsigned>>)
                return ItemOccurs<Char, IsOption, std::greater_equal<unsigned>>(this->m_name, this->m_count);
            else if constexpr (std::is_same_v<Comp, std::equal_to<unsigned>>)
                return ItemOccurs<Char, IsOption, std::not_equal_to<unsigned>>(this->m_name, this->m_count);
            else if constexpr (std::is_same_v<Comp, std::not_equal_to<unsigned>>)
                return ItemOccurs<Char, IsOption, std::equal_to<unsigned>>(this->m_name, this->m_count);
        }

        friend auto describe(const ItemOccurs & val)  {

            using Messages = Messages<Char>;
            
            auto typeName = IsOption ? Messages::option() : Messages::positionalArg();
            if constexpr (std::is_same_v<Comp, std::greater_equal<unsigned>>) {
                return format(Messages::itemMustBePresentGE(val.m_count), typeName, val.m_name, val.m_count);
            } else if constexpr (std::is_same_v<Comp, std::less_equal<unsigned>>) {
                return format(Messages::itemMustBePresentLE(val.m_count), typeName, val.m_name, val.m_count);
            } else if constexpr (std::is_same_v<Comp, std::greater<unsigned>>) {
                return format(Messages::itemMustBePresentG(val.m_count), typeName, val.m_name, val.m_count);
            } else if constexpr (std::is_same_v<Comp, std::less<unsigned>>) {
                return format(Messages::itemMustBePresentL(val.m_count), typeName, val.m_name, val.m_count);
            } else if constexpr (std::is_same_v<Comp, std::equal_to<unsigned>>) {
                return format(Messages::itemMustBePresentEQ(val.m_count), typeName, val.m_name, val.m_count);
            } else if constexpr (std::is_same_v<Comp, std::not_equal_to<unsigned>>) {
                return format(Messages::itemMustBePresentNEQ(val.m_count), typeName, val.m_name, val.m_count);
            } 
        }
    private:
        std::basic_string<Char> m_name;
        unsigned m_count;
    };


    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionPresent(S name) {
        return ItemOccurs<CharTypeOf<S>, true, std::greater<unsigned>>(name, 0);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionAbsent(S name) {
        return ItemOccurs<CharTypeOf<S>, true, std::equal_to<unsigned>>(name, 0);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionOccursAtLeast(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, true, std::greater_equal<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionOccursAtMost(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, true, std::less_equal<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionOccursMoreThan(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, true, std::greater<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionOccursLessThan(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, true, std::less<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionOccursExactly(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, true, std::equal_to<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto optionDoesntOccurExactly(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, true, std::not_equal_to<unsigned>>(name, count);
    }

    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalPresent(S name) {
        return ItemOccurs<CharTypeOf<S>, false, std::greater<unsigned>>(name, 0);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalAbsent(S name) {
        return ItemOccurs<CharTypeOf<S>, false, std::equal_to<unsigned>>(name, 0);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalOccursAtLeast(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, false, std::greater_equal<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalOccursAtMost(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, false, std::less_equal<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalOccursMoreThan(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, false, std::greater<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalOccursLessThan(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, false, std::less<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalOccursExactly(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, false, std::equal_to<unsigned>>(name, count);
    }
    ARGUM_MOD_EXPORTED template<StringLike S>
    auto positionalDoesntOccurExactly(S name, unsigned count) {
        return ItemOccurs<CharTypeOf<S>, false, std::not_equal_to<unsigned>>(name, count);
    }
    
}

#endif 
#ifndef HEADER_ARGUM_TOKENIZER_H_INCLUDED
#define HEADER_ARGUM_TOKENIZER_H_INCLUDED

#ifndef HEADER_ARGUM_EXPECTED_H_INCLUDED
#define HEADER_ARGUM_EXPECTED_H_INCLUDED

#ifndef HEADER_ARGUM_DATA_H_INCLUDED
#define HEADER_ARGUM_DATA_H_INCLUDED



namespace Argum {

    ARGUM_MOD_EXPORTED
    template<class Char>
    class BasicOptionNames final {
    
    public:
        using CharType = Char;
        using StringType = std::basic_string<Char>;
        using StringViewType = std::basic_string_view<Char>;
        using CharConstants = Argum::CharConstants<Char>;
    public:
        template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
        BasicOptionNames(First && first, Rest && ...rest):
            m_values{makeString(std::forward<First>(first)), makeString(std::forward<Rest>(rest))...} {
        }

        auto main() const -> const StringType & {
            return this->m_values[0];
        }

        auto all() const -> const std::vector<StringType> & {
            return this->m_values;
        }


    private:
        std::vector<StringType> m_values;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(OptionNames)

    enum class OptionArgumentKind  {
        None,
        Optional,
        Required
    };

    ARGUM_MOD_EXPORTED class Quantifier {
    public:
        static constexpr unsigned infinity = std::numeric_limits<unsigned>::max();

        constexpr Quantifier(unsigned val): m_min(val), m_max(val) {
        }
        constexpr Quantifier(unsigned min, unsigned max): m_min(min), m_max(max) {
            if (min > max)
                ARGUM_INVALID_ARGUMENT("min must be less or equal to max");
        }
        
        constexpr auto min() const {
            return this->m_min;
        }

        constexpr auto max() const {
            return this->m_max;
        }

    private:
        unsigned m_min = 0;
        unsigned m_max = 0;
    };

    ARGUM_MOD_EXPORTED inline constexpr Quantifier zeroOrOneTime  (0, 1);
    ARGUM_MOD_EXPORTED inline constexpr Quantifier neverOrOnce    = zeroOrOneTime;
    ARGUM_MOD_EXPORTED inline constexpr Quantifier oneTime        (1, 1);
    ARGUM_MOD_EXPORTED inline constexpr Quantifier once           = oneTime;
    ARGUM_MOD_EXPORTED inline constexpr Quantifier zeroOrMoreTimes(0, Quantifier::infinity);
    ARGUM_MOD_EXPORTED inline constexpr Quantifier oneOrMoreTimes (1, Quantifier::infinity);
    ARGUM_MOD_EXPORTED inline constexpr Quantifier onceOrMore     = oneOrMoreTimes;

    ARGUM_MOD_EXPORTED enum class Error {
        UnrecognizedOption = 1,
        AmbiguousOption,
        MissingOptionArgument,
        ExtraOptionArgument,
        ExtraPositional,
        ValidationError,
        ResponseFileError,

        Last = ResponseFileError,
        UserError = int(Last) + 100
    };

    ARGUM_MOD_EXPORTED
    template<class Char>
    class BasicParsingException : public std::exception {
    private:
        using StringsType = std::conditional_t<std::is_same_v<Char, char>, 
                                    std::tuple<std::string>,
                                    std::tuple<std::string, std::basic_string<Char>>>;
    public:
        auto message() const noexcept -> std::basic_string_view<Char> {
            return std::get<std::tuple_size_v<StringsType> - 1>(this->m_strings);
        }

        auto code() const noexcept -> Error {
            return m_code;
        }

        template<class Derived>
        requires(std::is_base_of_v<BasicParsingException, Derived>)
        auto as() const noexcept -> const Derived * {
            if (this->m_code == Derived::ErrorCode)
                return static_cast<const Derived *>(this);
            return nullptr;
        }

        template<class Derived>
        requires(std::is_base_of_v<BasicParsingException, Derived>)
        auto as() noexcept -> Derived * {
            return const_cast<Derived *>(const_cast<const BasicParsingException *>(this)->as<Derived>());
        }

        auto what() const noexcept -> const char * override {
            return std::get<0>(this->m_strings).c_str();
        }

        virtual auto clone() const & -> std::shared_ptr<BasicParsingException> = 0;
        virtual auto clone() & -> std::shared_ptr<BasicParsingException> = 0;
        virtual auto clone() && -> std::shared_ptr<BasicParsingException> = 0;
        [[noreturn]] virtual auto raise() const -> void = 0;
        
    protected:
        BasicParsingException(Error code, std::basic_string<Char> message) : 
            m_strings(BasicParsingException::initStrings(std::move(message))),
            m_code(code) {
        }
    private:
        static auto initStrings(std::basic_string<Char> && message) -> StringsType {
            if constexpr (std::is_same_v<Char, char>) {
                return StringsType(std::move(message));
            } else {
                return StringsType(toString<char>(message), std::move(message));
            }
        }
    private:
        StringsType m_strings;
        Error m_code;
    };
    
    ARGUM_DECLARE_FRIENDLY_NAMES(ParsingException)

    #define ARGUM_IMPLEMENT_EXCEPTION(type, base, code) \
        static constexpr Error ErrorCode = code; \
        auto clone() const & -> std::shared_ptr<base> override { \
            return std::make_shared<type>(*this); \
        } \
        auto clone() & -> std::shared_ptr<base> override { \
            return std::make_shared<type>(*this); \
        } \
        auto clone() && -> std::shared_ptr<base> override { \
            return std::make_shared<type>(std::move(*this)); \
        } \
        [[noreturn]] auto raise() const -> void override { ARGUM_RAISE_EXCEPTION(*this); }

}


#endif


namespace Argum {

    ARGUM_MOD_EXPORTED template<class T> using FailureType = std::in_place_type_t<T>;

    ARGUM_MOD_EXPORTED template<class T> inline constexpr FailureType<T> Failure{};

    ARGUM_MOD_EXPORTED
    template<class Char, class T>
    class [[nodiscard]] BasicExpected {
        template<class OtherChar, class Other> friend class BasicExpected;
    private:
        using ValueType = std::conditional_t<std::is_same_v<T, void>, std::monostate, T>;
    public:
        using ParsingException = BasicParsingException<Char>;
        using ParsingExceptionPtr = std::shared_ptr<ParsingException>;

        using ConstLValueReference = std::add_lvalue_reference_t<const T>;
        using LValueReference = std::add_lvalue_reference_t<T>;
        using RValueReference = std::add_rvalue_reference_t<T>;
    private:
        using ImplType = std::variant<ValueType, ParsingExceptionPtr>;
    public:
        BasicExpected() = default;

        template<class... Args>
        requires(std::is_constructible_v<T, Args &&...> && !std::is_same_v<T, void>)
        BasicExpected(Args && ...args): 
            m_impl(std::in_place_type_t<T>(), std::forward<Args>(args)...) {
        }

        BasicExpected(ParsingExceptionPtr err): m_impl(BasicExpected::validate(err)) {
        }

        template<class Exception, class... Args>
        requires(std::is_base_of_v<ParsingException, Exception>)
        BasicExpected(FailureType<Exception>, Args && ...args): 
            m_impl(std::make_shared<Exception>(std::forward<Args>(args)...)) {
        }

        template<class OtherT>
        requires(std::is_constructible_v<T, OtherT> || std::is_same_v<T, void>)
        BasicExpected(const BasicExpected<Char, OtherT> & other): 
            m_impl(std::visit([&](const auto & val) {
                    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(val)>, ParsingExceptionPtr>) {
                        return ImplType(val); 
                    } else if constexpr (!std::is_same_v<T, void>) {
                        return ImplType(val); 
                    } else {
                        return ImplType();
                    }
                }, other.m_impl)) {
        }

        auto operator*() const & -> ConstLValueReference {
            if constexpr (!std::is_same_v<T, void>)
                return *std::get_if<T>(&this->m_impl);
        }
        auto operator*() & -> LValueReference {
            if constexpr (!std::is_same_v<T, void>)
                return *std::get_if<T>(&this->m_impl);
        }
        auto operator*() && -> RValueReference {
            if constexpr (!std::is_same_v<T, void>)
                return std::move(*std::get_if<T>(&this->m_impl));
        }
        auto value() const & -> ConstLValueReference {
            return std::visit([](const auto & val) -> ConstLValueReference {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(val)>, ValueType>) {
                    if constexpr (!std::is_same_v<T, void>)
                        return val;
                } else {
                    BasicExpected::raise(val);
                }
            }, this->m_impl);
        }
        auto value() & -> LValueReference {
            return std::visit([](auto & val) -> LValueReference {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(val)>, ValueType>) {
                    if constexpr (!std::is_same_v<T, void>)
                        return val;
                } else {
                    BasicExpected::raise(val);
                }
            }, this->m_impl);
        }
        auto value() && -> RValueReference {
            return std::visit([](auto && val) -> RValueReference {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(val)>, ValueType>) {
                    if constexpr (!std::is_same_v<T, void>)
                        return std::forward<std::remove_cvref_t<decltype(val)>>(val);
                } else {
                    BasicExpected::raise(val);
                }
            }, std::move(this->m_impl));
        }

        template<class X = T>
        requires(!std::is_same_v<X, void>)
        auto operator->() const -> const T * {
            return std::get_if<T>(&this->m_impl);
        }
        template<class X = T>
        requires(!std::is_same_v<X, void>)
        auto operator->() -> T * {
            return std::get_if<T>(&this->m_impl);
        }

        auto error() const -> ParsingExceptionPtr {
            return std::visit([](const auto & val) -> ParsingExceptionPtr {
                if constexpr (std::is_same_v<std::remove_cvref_t<decltype(val)>, ParsingExceptionPtr>) {
                    return val;
                } else {
                    return ParsingExceptionPtr();
                }
            }, this->m_impl);
        }

        explicit operator bool() const {
            return std::holds_alternative<ValueType>(this->m_impl);
        }
        auto operator!() const -> bool {
            return !bool(*this);
        }
    private:
        static auto validate(ParsingExceptionPtr & ptr) -> ParsingExceptionPtr && {
            if (!ptr)
                ARGUM_INVALID_ARGUMENT("error must be non-null");
            return std::move(ptr);
        }
        [[noreturn]] static auto raise(const ParsingExceptionPtr & ptr) {
            if (ptr)
                ptr->raise();
            ARGUM_ALWAYS_ASSERT(!"accessing value of moved-out BasicExpected");
            abort();
        }
    private:
        ImplType m_impl;
    };

    ARGUM_MOD_EXPORTED template<class T> using Expected = BasicExpected<char, T>;
    ARGUM_MOD_EXPORTED template<class T> using WExpected = BasicExpected<wchar_t, T>;

    template<class T, class Char, class Ret, class... Args>
    concept IsHandler = std::is_invocable_v<std::remove_cvref_t<T>, Args...> && (
                            std::is_same_v<Ret, decltype(std::declval<std::remove_cvref_t<T>>()(std::declval<Args>()...))> ||
                            std::is_same_v<BasicExpected<Char, Ret>, decltype(std::declval<std::remove_cvref_t<T>>()(std::declval<Args>()...))>);

    template<class Char, class H, class Ret, class... Args>
    decltype(auto) adaptHandler(H && h) 
    requires(IsHandler<decltype(h), Char, Ret, Args...>) {

        using InHandlerType = std::remove_cvref_t<decltype(h)>;

        #ifdef ARGUM_USE_EXPECTED

            if constexpr (std::is_invocable_r_v<BasicExpected<Char, Ret>, InHandlerType, Args...>) {

                return std::forward<H>(h);

            } else {

                InHandlerType handler(std::forward<H>(h));

                return [=](Args ...args) -> BasicExpected<Char, Ret> {
                    if constexpr (std::is_same_v<Ret, void>) {
                        handler(args...);
                        return {};
                    } else {
                        return handler(args...);
                    }
                };
            }

        #else
            if constexpr (std::is_invocable_r_v<Ret, InHandlerType, Args...>) {

                return std::forward<H>(h);

            } else {

                InHandlerType handler(std::forward<H>(h));
                return [=](Args ...args) -> Ret {
                    return handler(args...).value();
                };
            }
        #endif
    }

    #ifdef ARGUM_USE_EXPECTED
        #define ARGUM_EXPECTED(c, type) BasicExpected<c, type>
        #define ARGUM_PROPAGATE_ERROR(expr) if (auto err = (expr).error()) { return err; }
        #define ARGUM_CHECK_RESULT_IMPL(temp, var, expr) decltype(auto) temp = (expr); ARGUM_PROPAGATE_ERROR(temp); var = *temp
        #define ARGUM_CHECK_RESULT(var, expr)  ARGUM_CHECK_RESULT_IMPL(ARGUM_UNIQUE_NAME(argum_check_result), var, expr)
        #define ARGUM_THROW(type, ...) return {Failure<type> __VA_OPT__(,) __VA_ARGS__}
        #define ARGUM_VOID_SUCCESS {}
    #else
        #define ARGUM_EXPECTED(c, x) x
        #define ARGUM_PROPAGATE_ERROR(expr) do { expr; } while(false)
        #define ARGUM_CHECK_RESULT(var, expr)  var = expr
        #define ARGUM_THROW(type, ...) throw type(__VA_ARGS__)
        #define ARGUM_VOID_SUCCESS
    #endif
}



#endif



namespace Argum {

    ARGUM_MOD_EXPORTED
    template<class Char>
    class BasicTokenizer final {

    private:
        using CharConstants = Argum::CharConstants<Char>;

        enum PrefixType {
            NotPrefix    = 0,
            ShortPrefix  = 1,
            LongPrefix   = 2,
            OptionStop   = 4
        };
        friend auto operator&(PrefixType lhs, PrefixType rhs) -> PrefixType {
            return PrefixType(std::underlying_type_t<PrefixType>(lhs) & std::underlying_type_t<PrefixType>(rhs));
        }
        friend auto operator&=(PrefixType & lhs, PrefixType rhs) -> PrefixType & {
            lhs = (lhs & rhs);
            return lhs;
        }
        friend auto operator|(PrefixType lhs, PrefixType rhs) -> PrefixType {
            return PrefixType(std::underlying_type_t<PrefixType>(lhs) | std::underlying_type_t<PrefixType>(rhs));
        }
        friend auto operator|=(PrefixType & lhs, PrefixType rhs) -> PrefixType & {
            lhs = (lhs | rhs);
            return lhs;
        }
        using PrefixId = unsigned;
        using NameIndex = unsigned;

    public:
        using CharType = Char;
        using StringType = std::basic_string<Char>;
        using StringViewType = std::basic_string_view<Char>;
        using OptionNames = BasicOptionNames<Char>;

        struct OptionToken {
            unsigned argIdx;                        //index of the argument containing the token in command line
            unsigned idx;                           //index of the option
            StringType usedName;                    //specific option name used
            std::optional<StringViewType> argument; //argument if included with the otpion via = syntax

            auto operator==(const OptionToken & rhs) const -> bool = default;
            auto operator!=(const OptionToken & rhs) const -> bool = default;
        };

        struct ArgumentToken {
            unsigned argIdx;                        //index of the argument containing the token in command line
            StringViewType value;                   //text of the argument

            auto operator==(const ArgumentToken & rhs) const -> bool = default;
            auto operator!=(const ArgumentToken & rhs) const -> bool = default;
        };

        struct OptionStopToken {
            unsigned argIdx;                        //index of the argument containing the token in command line

            auto operator==(const OptionStopToken & rhs) const -> bool = default;
            auto operator!=(const OptionStopToken & rhs) const -> bool = default;
        };

        struct UnknownOptionToken {
            unsigned argIdx;                        //index of the argument containing the token in command line
            StringType name;                        //option name
            std::optional<StringViewType> argument; //argument if included with the otpion via = syntax

            auto operator==(const UnknownOptionToken & rhs) const -> bool = default;
            auto operator!=(const UnknownOptionToken & rhs) const -> bool = default;
        };

        struct AmbiguousOptionToken {
            unsigned argIdx;                        //index of the argument containing the token in command line
            StringType name;                        //option name
            std::optional<StringViewType> argument; //argument if included with the otpion via = syntax
            std::vector<StringType> possibilities;  //possible completions

            auto operator==(const AmbiguousOptionToken & rhs) const -> bool = default;
            auto operator!=(const AmbiguousOptionToken & rhs) const -> bool = default;
        };


        enum TokenResult {
            Continue   = 0,
            StopAfter  = 0b10,
            StopBefore = 0b11
        };

        class Settings {
            friend BasicTokenizer;
        public:
            template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
            auto addLongPrefix(First && first, Rest && ...rest) & -> Settings & {
                std::initializer_list<StringType> values = {makeString(std::forward<First>(first)), makeString(std::forward<Rest>(rest))...};
                for(auto & value: values) {
                    auto [it, inserted] = m_prefixes.add(std::move(value), m_lastPrefixId);
                    auto & type = m_prefixTypes[it->value()];
                    if (!inserted) {
                        if ((type & ShortPrefix) == ShortPrefix)
                            ARGUM_INVALID_ARGUMENT("the same prefix cannot be used for long and short options");
                    }
                    type |= LongPrefix;
                }
                ++m_lastPrefixId;
                return *this;
            }

            template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
            auto addLongPrefix(First && first, Rest && ...rest) && -> Settings && {
                return std::move(static_cast<Settings *>(this)->addLongPrefix(std::forward<First>(first), std::forward<Rest>(rest)...));
            }

            template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
            auto addShortPrefix(First && first, Rest && ...rest) & -> Settings & {
                std::initializer_list<StringType> values = {makeString(std::forward<First>(first)), makeString(std::forward<Rest>(rest))...};
                for(auto & value: values) {
                    auto [it, inserted] = m_prefixes.add(std::move(value), m_lastPrefixId);
                    auto & type = m_prefixTypes[it->value()];
                    if (!inserted) {
                        //the same prefix cannot be used for long and short
                        if ((type & LongPrefix) == LongPrefix) 
                            ARGUM_INVALID_ARGUMENT("the same prefix cannot be used for long and short options");
                    }
                    type |= ShortPrefix;
                }
                ++m_lastPrefixId;
                return *this;
            }

            template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
            auto addShortPrefix(First && first, Rest && ...rest) && -> Settings && {
                return std::move(static_cast<Settings *>(this)->addShortPrefix(std::forward<First>(first), std::forward<Rest>(rest)...));
            }

            template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
            auto addOptionTerminator(First && first, Rest && ...rest) & -> Settings & {
                std::initializer_list<StringType> values = {makeString(std::forward<First>(first)), makeString(std::forward<Rest>(rest))...};
                for(auto & value: values) {
                    auto [it, inserted] = m_prefixes.add(std::move(value), m_lastPrefixId);
                    auto & type = m_prefixTypes[it->value()];
                    type |= OptionStop;
                }
                ++m_lastPrefixId;
                return *this;
            }

            template<StringConvertibleOf<CharType> First, StringConvertibleOf<CharType>... Rest>
            auto addOptionTerminator(First && first, Rest && ...rest) && -> Settings && {
                return std::move(static_cast<Settings *>(this)->addOptionTerminator(std::forward<First>(first), std::forward<Rest>(rest)...));
            }

            auto addValueDelimiter(CharType c) & -> Settings & {
                
                auto [it, inserted] = m_valueDelimiters.insert(c);
                if (!inserted)
                    ARGUM_INVALID_ARGUMENT("duplicate delimiter");
                return *this;
            }

            auto addValueDelimiter(CharType c) && -> Settings && {
                return std::move(static_cast<Settings *>(this)->addValueDelimiter(c));
            }

            auto allowAbbreviation(bool value) & -> Settings & {
                m_allowAbrreviation = value;
                return *this;
            }

            auto allowAbbreviation(bool value) && -> Settings && {
                return std::move(static_cast<Settings *>(this)->allowAbbreviation(value));
            }

            static auto commonUnix() -> Settings {
                Settings ret;
                ret.addLongPrefix(CharConstants::doubleDash)
                   .addShortPrefix(CharConstants::dash)
                   .addOptionTerminator(CharConstants::doubleDash)
                   .addValueDelimiter(CharConstants::assignment);
                return ret;
            }

            static auto unixLongOnly() -> Settings {
                Settings ret;
                ret.addLongPrefix(CharConstants::doubleDash, CharConstants::dash)
                   .addOptionTerminator(CharConstants::doubleDash)
                   .addValueDelimiter(CharConstants::assignment);
                return ret;
            }

            static auto windowsShort() -> Settings {
                Settings ret;
                ret.addShortPrefix(CharConstants::slash, CharConstants::dash)
                   .addOptionTerminator(CharConstants::doubleDash)
                   .addValueDelimiter(CharConstants::colon);
                return ret;
            }

            static auto windowsLong() -> Settings {
                Settings ret;
                ret.addLongPrefix(CharConstants::slash, CharConstants::dash, CharConstants::doubleDash)
                   .addOptionTerminator(CharConstants::doubleDash)
                   .addValueDelimiter(CharConstants::colon);
                return ret;
            }
        private: 
            FlatMap<StringType, PrefixId> m_prefixes;
            FlatMap<PrefixId, PrefixType> m_prefixTypes;
            FlatSet<CharType> m_valueDelimiters;
            PrefixId m_lastPrefixId = 0;
            bool m_allowAbrreviation = true;
        };
        
    public:
        BasicTokenizer():
            BasicTokenizer(Settings::commonUnix()) {
        }

        BasicTokenizer(Settings settings) {

            this->m_prefixes = std::move(settings.m_prefixes);
            this->m_prefixTypes = std::move(settings.m_prefixTypes);
            this->m_valueDelimiters = std::move(settings.m_valueDelimiters);
            this->m_allowAbrreviation = settings.m_allowAbrreviation;
        }

        auto add(const OptionNames & names)  {

            auto currentIndex = unsigned(this->m_names.size());
            for(auto & opt: names.all()) {
                auto findResult = this->findLongestPrefix(opt);
                if (!findResult)
                    ARGUM_INVALID_ARGUMENT("option must start with a valid prefix");
                if (findResult->size == opt.size())
                    ARGUM_INVALID_ARGUMENT("option must have more than a prefix");

                if ((findResult->type & LongPrefix) == LongPrefix) {
                    this->add(this->m_longs[findResult->index], opt.substr(findResult->size), currentIndex);
                } else if ((findResult->type & ShortPrefix) == ShortPrefix) {
                    if (opt.size() - findResult->size == 1)
                        this->add(this->m_singleShorts[findResult->index], opt[findResult->size], currentIndex);
                    else
                        this->add(this->m_multiShorts[findResult->index], opt.substr(findResult->size), currentIndex);
                } else {
                    ARGUM_INVALID_ARGUMENT("option is neither short nor long with currently defined prefixes");
                }
            }
            this->m_names.emplace_back(names.main());
        }

        template<ArgIterator<CharType> It, class Func>
        auto tokenize(It argFirst, It argLast, Func && handler) const -> ARGUM_EXPECTED(CharType, std::vector<StringType>)  {

            static_assert(IsHandler<Func, CharType, TokenResult, OptionToken>, "handler must handle OptionToken and return correct return type");
            static_assert(IsHandler<Func, CharType, TokenResult, ArgumentToken>, "handler must handle ArgumentToken and return correct return type");
            static_assert(IsHandler<Func, CharType, TokenResult, OptionStopToken>, "handler must handle OptionStopToken and return correct return type");
            static_assert(IsHandler<Func, CharType, TokenResult, UnknownOptionToken>, "handler must handle UnknownOptionToken and return correct return type");
            static_assert(IsHandler<Func, CharType, TokenResult, AmbiguousOptionToken>, "handler must handle AmbiguousOptionToken and return correct return type");

            bool noMoreOptions = false;
            std::vector<StringType> rest;
            StringViewType arg; //last argument being processed
            unsigned consumed = 0; //how much of the arg was consumed before invoking a handler
            unsigned unconsumedPrefixSize = 0; //size of the arg option prefix if it is an option

            for(unsigned argIdx = 0; argFirst != argLast; ++argFirst, ++argIdx) {
                arg = *argFirst;
                consumed = 0;
                unconsumedPrefixSize = 0;

                std::optional<TokenResult> result;
                if (!noMoreOptions) {

                    if (auto prefixFindResult = this->findLongestPrefix(arg)) {
                    
                        auto type = prefixFindResult->type;
                    
                        if (prefixFindResult->size == arg.size()) {
                            
                            if ((type & OptionStop) == OptionStop) {
                                noMoreOptions = true;
                                ARGUM_CHECK_RESULT(result, this->callHandler(std::forward<Func>(handler), OptionStopToken{argIdx}));
                                if (result == TokenResult::StopAfter)
                                    consumed = unsigned(arg.size());
                            }
                            
                        } else {
                            if ((type & LongPrefix) == LongPrefix) {
                                ARGUM_CHECK_RESULT(result, this->handleLongPrefix(argIdx, arg,
                                                                         prefixFindResult->index, prefixFindResult->size,
                                                                         std::forward<Func>(handler)));
                                if (result == TokenResult::StopAfter)
                                    consumed = unsigned(arg.size());
                            } else if ((type & ShortPrefix) == ShortPrefix) {
                                ARGUM_CHECK_RESULT(result, this->handleShortPrefix(argIdx, arg,
                                                                          prefixFindResult->index, prefixFindResult->size,
                                                                          consumed, handler));
                                unconsumedPrefixSize = prefixFindResult->size;
                            }
                        }
                    }
                }
                
                if (!result) {
                    ARGUM_CHECK_RESULT(result, this->callHandler(std::forward<Func>(handler), ArgumentToken{argIdx, arg}));
                    if (result == TokenResult::StopAfter)
                        consumed = unsigned(arg.size());
                }

                if (*result != TokenResult::Continue)
                    break;
            }

            if (argFirst != argLast) {
                
                if (consumed == arg.size()) {
                    ++argFirst;
                } else if (consumed != 0) {
                    rest.emplace_back(StringType(arg.substr(0, unconsumedPrefixSize)) + StringType(arg.substr(consumed)));
                    ++argFirst;
                }

                for ( ; argFirst != argLast; ++argFirst) {
                    rest.emplace_back(*argFirst);
                }
            }
            return rest;
        }

        auto isOptionNameLong(const StringViewType name) const -> bool {
            
            auto findResult = this->findLongestPrefix(name);
            if (!findResult)
                ARGUM_INVALID_ARGUMENT("option must start with a valid prefix");
            if (findResult->size == name.size())
                ARGUM_INVALID_ARGUMENT("option must have more than a prefix");

            if ((findResult->type & LongPrefix) == LongPrefix) {
                return true;
            } else if ((findResult->type & ShortPrefix) == ShortPrefix) {
                return false;
            } 
                
            ARGUM_INVALID_ARGUMENT("option is neither short nor long with currently defined prefixes");
        }

    private:
        template<class Token, class Handler>
        auto callHandler(Handler && handler, Token && token) const -> ARGUM_EXPECTED(CharType, TokenResult) {

            return adaptHandler<CharType, 
                                Handler, 
                                TokenResult,
                                std::remove_cvref_t<Token>>(std::forward<Handler>(handler))(std::forward<Token>(token));
        }


        template<class Func>
        auto handleLongPrefix(unsigned argIdx, 
                              StringViewType option, 
                              PrefixId prefixId,
                              unsigned nameStart,
                              Func && handler) const -> ARGUM_EXPECTED(CharType, TokenResult) {
            
            auto [name, arg] = this->splitDelimitedArgument(option, nameStart);
            if (name.size() == 0)
                return this->callHandler(std::forward<Func>(handler), ArgumentToken{argIdx, option});

            StringType usedName(option.begin(), 
#if defined(_MSC_VER) && _ITERATOR_DEBUG_LEVEL > 0
                option.begin() + nameStart + name.size()
#else
                name.end()
#endif
            );

            auto mapIt = this->m_longs.find(prefixId);
            if (mapIt == this->m_longs.end()) {
                return this->callHandler(std::forward<Func>(handler), UnknownOptionToken{argIdx, std::move(usedName), std::move(arg)});
            }
            auto & longsMap = mapIt->value();

            if (this->m_allowAbrreviation) {
                const auto & [first, last] = findMatchOrMatchingPrefixRange(longsMap, name);
                if (last - first == 1) {
                    return this->callHandler(std::forward<Func>(handler), OptionToken{argIdx, first->value(), std::move(usedName), std::move(arg)});
                } else if (last != first) {
                    StringType actualPrefix(option.substr(0, nameStart));
                    std::vector<StringType> candidates(last - first);
                    std::transform(first, last, candidates.begin(), [&](const auto & p) {return actualPrefix + p.key(); });
                    return this->callHandler(std::forward<Func>(handler), AmbiguousOptionToken{argIdx, std::move(usedName), std::move(arg), candidates});
                }
            } else {
                auto it = longsMap.find(name);
                if (it != longsMap.end())
                    return this->callHandler(std::forward<Func>(handler), OptionToken{argIdx, it->value(), std::move(usedName), std::move(arg)});
            }
            
            if (auto maybeToken = this->matchNumber(argIdx, option, nameStart)) {
                return this->callHandler(std::forward<Func>(handler), *maybeToken);
            }
            return this->callHandler(std::forward<Func>(handler), UnknownOptionToken{argIdx, std::move(usedName), std::move(arg)});
        }

        template<class Func>
        auto handleShortPrefix(unsigned argIdx, 
                               StringViewType option, 
                               PrefixId prefixId,
                               unsigned nameStart,
                               unsigned & consumed,
                               Func && handler) const -> ARGUM_EXPECTED(CharType, TokenResult) {

            if (auto maybeResult = this->handleShortOption(argIdx, option, prefixId, nameStart, consumed, handler)) {
                return std::move(*maybeResult);
            } 

            TokenResult result;
            if (auto maybeToken = this->matchNumber(argIdx, option, nameStart)) {
                ARGUM_CHECK_RESULT(result, this->callHandler(std::forward<Func>(handler), *maybeToken));
            } else {
                ARGUM_CHECK_RESULT(result, this->callHandler(std::forward<Func>(handler), UnknownOptionToken{argIdx, StringType(option), std::nullopt}));
            }

            if (result == TokenResult::StopAfter)
                consumed = unsigned(option.size());
            return result;
        }

        template<class Func>
        auto handleShortOption(unsigned argIdx, 
                               StringViewType option, 
                               PrefixId prefixId,
                               unsigned nameStart, 
                               unsigned & consumed,
                               Func && handler) const -> std::optional<ARGUM_EXPECTED(CharType, TokenResult)> {

            StringViewType chars = option.substr(nameStart);
            assert(!chars.empty());

            std::optional<unsigned> singleLetterNameIdx;
            
            auto mapIt = this->m_singleShorts.find(prefixId);
            if (mapIt != this->m_singleShorts.end()) {
                auto it = mapIt->value().find(chars[0]);
                if (it != mapIt->value().end())
                    singleLetterNameIdx = it->value();
            }

            if (chars.size() > 1 || !singleLetterNameIdx) {

                ARGUM_CHECK_RESULT(auto maybeResult, 
                    this->handleMultiShortOption(argIdx, option, prefixId, nameStart, 
                                                 singleLetterNameIdx.has_value(), std::forward<Func>(handler)));
                if (maybeResult) {
                    if (*maybeResult == TokenResult::StopAfter)
                        consumed = unsigned(option.size());
                    return *maybeResult;
                }
            }

            if (!singleLetterNameIdx)
                return std::nullopt;
            
            auto & shortsMap = mapIt->value();
            consumed = nameStart;
            auto actualPrefix = StringType(option.substr(0, nameStart));
            do {
                
                auto currentIdx = *singleLetterNameIdx;
                auto usedName = actualPrefix + chars[0];
                std::optional<StringViewType> arg;

                unsigned charsConsumed = 1;
                if (chars.size() > 1) {

                    auto it = shortsMap.find(chars[1]);
                    if (it == shortsMap.end()) {
                        arg = chars.substr(1);
                        charsConsumed = unsigned(chars.size());
                    }
                    else {
                        singleLetterNameIdx = it->value();
                    }
                }

                ARGUM_CHECK_RESULT(auto result, 
                    this->callHandler(std::forward<Func>(handler), 
                                      OptionToken{argIdx, currentIdx, std::move(usedName), std::move(arg)}));
                if (result != TokenResult::Continue) {
                    if (result == TokenResult::StopAfter)
                        consumed += charsConsumed;
                    return result;
                }

                chars.remove_prefix(charsConsumed);
                consumed += charsConsumed;
                
            } while(!chars.empty());

            return TokenResult::Continue;
        }

        template<class Func>
        auto handleMultiShortOption(unsigned argIdx, 
                                    StringViewType option,
                                    PrefixId prefixId,
                                    unsigned nameStart,
                                    bool mustMatchExact,
                                    Func && handler) const -> ARGUM_EXPECTED(CharType, std::optional<TokenResult>) {

            auto [name, arg] = this->splitDelimitedArgument(option, nameStart);
            if (name.size() == 0)
                return this->callHandler(std::forward<Func>(handler), ArgumentToken{argIdx, option});

            auto mapIt = this->m_multiShorts.find(prefixId);
            if (mapIt == this->m_multiShorts.end()) {
                return std::nullopt;
            }
            auto & multiShortsMap = mapIt->value();

            if (this->m_allowAbrreviation) {
                const auto & [first, last] = findMatchOrMatchingPrefixRange(multiShortsMap, name);
                if (last != first) {
                    StringType usedName(option.data(), name.data() + name.size());
                    if (last - first == 1) {
                        if (!mustMatchExact || first->key() == name) {
                            return this->callHandler(std::forward<Func>(handler), OptionToken{argIdx, first->value(), std::move(usedName), std::move(arg)});
                        } else {
                            std::vector<StringType> candidates = {
                                StringType(option.substr(0, nameStart + 1)),
                                StringType(option.substr(0, nameStart)) + first->key()
                            };
                            return this->callHandler(std::forward<Func>(handler), AmbiguousOptionToken{argIdx, std::move(usedName), std::move(arg), candidates});
                        }
                    } else  {
                        StringType actualPrefix(option.substr(0, nameStart));
                        std::vector<StringType> candidates;
                        if (mustMatchExact) {
                            candidates.reserve(1 + (last - first));
                            candidates.push_back(StringType(option.substr(0, nameStart + 1)));
                        } else {
                            candidates.reserve(last - first);
                        }
                        std::transform(first, last, std::back_inserter(candidates), [&](const auto & p) { return actualPrefix + p.key(); });
                        return this->callHandler(std::forward<Func>(handler), AmbiguousOptionToken{argIdx, std::move(usedName), std::move(arg), candidates});
                    }
                }
            } else {
                auto it = multiShortsMap.find(name);
                if (it != multiShortsMap.end()) {
                    StringType usedName(option.data(), name.data() + name.size());
                    return this->callHandler(std::forward<Func>(handler), OptionToken{argIdx, it->value(), std::move(usedName), std::move(arg)});
                }
            }
        

            return std::nullopt;
        }

        auto splitDelimitedArgument(StringViewType option, unsigned nameStart) const -> std::pair<StringViewType, std::optional<StringViewType>> {

            StringViewType name = option.substr(nameStart);
            std::optional<StringViewType> arg;

            for(auto symbol: this->m_valueDelimiters) {
                if (auto assignPos = option.find(symbol, nameStart); 
                    assignPos != option.npos && assignPos - nameStart < name.size()) {

                    name = option.substr(nameStart, assignPos - nameStart);
                    arg = option.substr(assignPos + 1);
                }
            }
            return {std::move(name), std::move(arg)};
        }

        auto matchNumber(unsigned argIdx, 
                         StringViewType option, 
                         unsigned /*nameStart*/) const -> std::optional<ArgumentToken> {

            //we assume that option is null terminated here!

            errno = 0;
            const CharType * p = &option[0];
            CharType * pEnd;
            CharConstants::toLongLong(p, &pEnd, 0);
            if (size_t(pEnd - p) == option.size() && errno != ERANGE)
                return ArgumentToken{argIdx, option};

            long double dres = CharConstants::toLongDouble(p, &pEnd);
            if (size_t(pEnd - p) == option.size() && dres != HUGE_VALL)
                return ArgumentToken{argIdx, option};
                
            return std::nullopt;
        }

        template<class Value, class Arg>
        static auto add(FlatMap<Value, unsigned> & map, Arg arg, unsigned idx) -> void {

            auto [it, inserted] = map.add(arg, idx);
            if (!inserted)
                ARGUM_INVALID_ARGUMENT("duplicate option");
        }

        struct PrefixFindResult {
            unsigned index; //index of the prefix in m_prefixes;
            unsigned size;  //size of specific matched value of the prefix
            PrefixType type;  
        };

        auto findLongestPrefix(StringViewType arg) const -> std::optional<PrefixFindResult> {

            std::optional<PrefixFindResult> ret;
            for(auto entry: this->m_prefixes) {
                auto & prefix = entry.key();
                auto prefixId = entry.value();
                
                auto typeIt = this->m_prefixTypes.find(prefixId);
                assert(typeIt != this->m_prefixTypes.end());

                if (matchPrefix(arg, StringViewType(prefix))) {
                    if (!ret || prefix.size() > ret->size)
                        ret = {prefixId, unsigned(prefix.size()), typeIt->value()};
                }
                
            }
            return ret;
        }
    private:
        FlatMap<StringType, PrefixId> m_prefixes;
        FlatMap<PrefixId, PrefixType> m_prefixTypes;
        FlatSet<CharType> m_valueDelimiters;
        std::vector<StringType> m_names;

        
        FlatMap<PrefixId, FlatMap<CharType, NameIndex>> m_singleShorts;
        FlatMap<PrefixId, FlatMap<StringType, NameIndex>> m_multiShorts;
        FlatMap<PrefixId, FlatMap<StringType, NameIndex>> m_longs;

        bool m_allowAbrreviation = true;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(Tokenizer)
}

#endif


#ifndef HEADER_ARGUM_PARTITIONER_H_INCLUDED
#define HEADER_ARGUM_PARTITIONER_H_INCLUDED




namespace Argum {

    // Given:
    //   1. A length of a sequence N and
    //   2. A list of M pairs {A,B} each denoting a range of minimum A and maximum B elements (A >=0 && A <= B) 
    // This class finds a partition of range {N1, N2, N3, ..., Nm+1} such as:
    //   a. N1 + N2 + N3 + ... + Nm + Nm+1 == N
    //   b. Ai <= Ni <= Bi for i <= m
    // The last partion: Nm+1 is the "remainder", if any.
    // The partitioning is "greedy" - each range consumes as much as it can from left to right
    template<class S>
    requires(std::is_integral_v<S>)
    class Partitioner {
    public:
        using SizeType = S;
        
        static constexpr auto infinity = std::numeric_limits<SizeType>::max();
    
    public:
        auto addRange(SizeType a, SizeType b) -> void {
            if constexpr (std::is_signed_v<SizeType>)
                assert(a >= 0);
            assert(a <= b);

            SizeType length;
            if (b != infinity)
                length = b - a;
            else
                length = infinity;
            this->m_ranges.push_back({a, length});
            this->m_minimumExpected += a;
        }

        //returns M + 1: e.g. the number of added ranges + 1
        auto paritionsCount() -> size_t {
            return this->m_ranges.size() + 1;
        }

        //The minimum size of sequence that can be partitioned
        auto minimumSequenceSize() const {
            return this->m_minimumExpected;
        }

        //Returns false if n < minimumSequenceSize()
        auto partition(SizeType n) -> std::optional<std::vector<SizeType>> {
            if (n < this->m_minimumExpected)
                return std::nullopt;
            //this "rebases" the sequence so now all we need is to match ranges of {0, length1}, {0, length2},...
            n -= this->m_minimumExpected;
            
            //because all ranges have minimum 0 and matching is greedy we can just assign
            //smaller of length and remaining n to each adding its minimum
            std::vector<SizeType> results(this->paritionsCount());
            std::transform(this->m_ranges.begin(), this->m_ranges.end(), results.begin(), [&n](const auto range) {
                auto length = std::min(n, range.second);
                n -= length;
                return SizeType(range.first + length);
            });
            results.back() = n;
            return results;
        }

    private:
        std::vector<std::pair<SizeType, SizeType>> m_ranges;
        SizeType m_minimumExpected = 0;
    };
}

#endif
#ifndef HEADER_ARGUM_COMMAND_LINE_H_INCLUDED
#define HEADER_ARGUM_COMMAND_LINE_H_INCLUDED

#ifndef HEADER_ARGUM_SIMPLE_FILE_H_INCLUDED
#define HEADER_ARGUM_SIMPLE_FILE_H_INCLUDED




namespace Argum {

    class SimpleFile {
    public:
        SimpleFile(const std::filesystem::path & path, const char * mode, std::error_code & ec) {
        #ifndef _MSC_VER
            m_fp = fopen(path.string().c_str(), mode);
        #else
            m_fp = _wfopen(path.native().c_str(), toString<wchar_t>(mode).c_str());
        #endif
            if (!m_fp) {
                int err = errno;
                ec = std::make_error_code(static_cast<std::errc>(err));
            }
        }
        SimpleFile(SimpleFile && src) {
            m_fp = src.m_fp;
            src.m_fp = nullptr;
        }
        auto operator=(SimpleFile && src) -> SimpleFile & {
            if (m_fp)
                fclose(m_fp);
            m_fp = nullptr;
            m_fp = src.m_fp;
            return *this;
        }
        ~SimpleFile() {
            if (m_fp)
                fclose(m_fp);
        }
        SimpleFile(const SimpleFile &) = delete;
        SimpleFile & operator=(const SimpleFile &) = delete;

        operator bool() const {
            return m_fp != nullptr;
        }
        auto eof() const -> bool {
            return feof(m_fp) != 0;
        }

        template<Character Char>
        std::basic_string<Char> readLine(std::error_code & ec) {
            std::basic_string<Char> buf;
            for ( ; ; ) {
                Char c;
                if constexpr (std::is_same_v<Char, char>) {
                    int res = fgetc(m_fp);
                    if (res == EOF) {
                        checkError(ec);
                        break;
                    }
                    c = Char(res);
                } else if constexpr (std::is_same_v<Char, wchar_t>) {
                    wint_t res = fgetwc(m_fp);
                    if (res == WEOF) {
                        checkError(ec);
                        break;
                    }
                    c = Char(res);
                }
                
                if (c == CharConstants<Char>::endl)
                    break;
                buf += c;
            }
            return buf;
        }
    private:
        auto checkError(std::error_code & ec) const -> bool{
            if (ferror(m_fp)) {
                int err = errno;
                 ec = std::make_error_code(static_cast<std::errc>(err));
                 return true;
            }
            return false;
        }
    private:
        FILE * m_fp;
    };
}

#endif

namespace Argum {

    ARGUM_MOD_EXPORTED 
    template<class Char>
    auto makeArgSpan(int argc, Char ** argv) -> std::span<const Char *> {
        if (argc > 0)
            return std::span(const_cast<const Char **>(argv + 1), size_t(argc - 1));
        return std::span(const_cast<const Char **>(argv), size_t(0));
    }

    ARGUM_MOD_EXPORTED
    template<class Char>
    class BasicResponseFileReader {
    public:
        using CharType = Char;
        using StringType = std::basic_string<Char>;
        using StringViewType = std::basic_string_view<Char>;

        struct Exception : public BasicParsingException<CharType> {
            ARGUM_IMPLEMENT_EXCEPTION(Exception, BasicParsingException<CharType>, Error::ResponseFileError)

            Exception(const std::filesystem::path & filename_, std::error_code error_): 
                BasicParsingException<CharType>(ErrorCode,
                                                format(Messages<CharType>::errorReadingResponseFile(), filename_.native(), error_.message())),
                filename(filename_),
                error(error_) {
            }
            std::filesystem::path filename;
            std::error_code error;
        };

    public:
        BasicResponseFileReader(char prefix) {
            m_prefixes.emplace_back(1, prefix);
        }
        BasicResponseFileReader(StringType prefix): m_prefixes(1, prefix) {
            m_prefixes.emplace_back(std::move(prefix));
        }
        BasicResponseFileReader(std::initializer_list<StringType> prefixes): m_prefixes(prefixes) {
        }

        auto expand(int argc, CharType ** argv) -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            return expand(makeArgSpan(argc, argv));
        }

        template<class Splitter>
        auto expand(int argc, CharType ** argv, Splitter && splitter) -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            return expand(makeArgSpan(argc, argv), std::forward<Splitter>(splitter));
        }

        template<ArgRange<CharType> Args>
        auto expand(const Args & args) -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            return expand(args, [](StringType && str, auto dest) {
                trimInPlace(str);
                if (str.empty())
                    return;
                *dest = std::move(str);
            });
        }

        template<ArgRange<Char> Args, class Splitter>
        auto expand(const Args & args, Splitter && splitter) -> ARGUM_EXPECTED(CharType, std::vector<StringType>)
            requires(std::is_invocable_v<decltype(splitter), StringType &&, std::back_insert_iterator<std::vector<StringType>>>) {
            
            std::vector<StringType> ret;
            std::stack<StackEntry> stack;

            for(StringViewType arg: args) {

                ARGUM_PROPAGATE_ERROR(this->handleArg(arg, ret, stack, splitter));

                while(!stack.empty()) {

                    auto & entry = stack.top();
                    if (entry.current == entry.items.end()) {
                        stack.pop();
                        continue;
                    }
                    ARGUM_PROPAGATE_ERROR(this->handleArg(*entry.current, ret, stack, std::forward<Splitter>(splitter)));
                    ++entry.current;
                }
            }
            return ret;
        }
    private:
        struct StackEntry {
            std::vector<StringType> items;
            typename std::vector<StringType>::const_iterator current;
        };

        template<class Splitter>
        auto handleArg(StringViewType arg, std::vector<StringType> & dest, 
                       std::stack<StackEntry> & stack, 
                       Splitter && splitter) -> ARGUM_EXPECTED(CharType, void) {

            auto foundIt = std::find_if(this->m_prefixes.begin(), this->m_prefixes.end(), [arg](const StringViewType & prefix) {
                return matchStrictPrefix(arg, prefix);
            });

            if (foundIt != this->m_prefixes.end()) {
                auto filename = arg.substr(foundIt->size());
                StackEntry nextEntry;
                ARGUM_PROPAGATE_ERROR(this->readResponseFile(filename, nextEntry.items, std::forward<Splitter>(splitter)));
                nextEntry.current = nextEntry.items.cbegin();
                stack.emplace(std::move(nextEntry));

            } else {
                dest.emplace_back(arg);
            }
            return ARGUM_VOID_SUCCESS;
        }
    
        template<class Splitter>
        static auto readResponseFile(StringViewType filename, std::vector<StringType> & dest, 
                                     Splitter && splitter) -> ARGUM_EXPECTED(CharType, void){

            std::filesystem::path path(filename);
            std::error_code error;
            SimpleFile file(path, "r", error);
            if (!file)
                ARGUM_THROW(Exception, path, error);

            do {
                StringType line = file.readLine<CharType>(error);
                if (error)
                    ARGUM_THROW(Exception, path, error);
                
                if (!line.empty())
                    splitter(std::move(line), std::back_inserter(dest));
            } while(!file.eof());

            return ARGUM_VOID_SUCCESS;
        }
    private:
        std::vector<StringType> m_prefixes;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(ResponseFileReader)

}

#endif 
#ifndef HEADER_ARGUM_HELP_FORMATTER_H_INCLUDED
#define HEADER_ARGUM_HELP_FORMATTER_H_INCLUDED

#ifndef HEADER_ARGUM_COLOR_H_INCLUDED
#define HEADER_ARGUM_COLOR_H_INCLUDED


namespace Argum {

    ARGUM_MOD_EXPORTED
    enum class Color : unsigned {
        normal = 0,
        bold = 1,
        faint = 2,
        italic = 3,
        underline = 4,
        reverse = 7,

        black = 30,
        red = 31,
        green = 32,
        yellow = 33,
        blue = 34,
        magenta = 35,
        cyan = 36,
        white = 37,

        bg_black = 40,
        bg_red = 41,
        bg_green = 42,
        bg_yellow = 43,
        bg_blue = 44,
        bg_magenta = 45,
        bg_cyan = 46,
        bg_white = 47,

        bright_black = 90,
        grey = bright_black,
        bright_red = 91,
        bright_green = 92,
        bright_yellow = 93,
        bright_blue = 94,
        bright_magenta = 95,
        bright_cyan = 96,
        bright_white = 97,

        bg_bright_black = 100,
        bg_grey = bg_bright_black,
        bg_bright_red = 101,
        bg_bright_green = 102,
        bg_bright_yellow = 103,
        bg_bright_blue = 104,
        bg_bright_magenta = 105,
        bg_bright_cyan = 106,
        bg_bright_white = 107
    };

    ARGUM_MOD_EXPORTED
    template<Character Char>
    class BasicColorMaker {
    private:
        using mytype = BasicColorMaker;

        template<size_t N>
        struct Str {
            Char wrapped[N];
        };

        template<size_t N>
        static constexpr auto append(Str<N> prev, Char c) {
            Str<N + 1> ret;
            std::copy(std::begin(prev.wrapped), std::begin(prev.wrapped) + N - 1, std::begin(ret.wrapped));
            ret.wrapped[N - 1] = c;
            ret.wrapped[N] = 0;
            return ret;
        }

        template<Color C, size_t N>
        static constexpr auto encode(Str<N> prev) {
            if constexpr (unsigned(C) < 10) {
                auto ret = mytype::append(prev, Char(CharConstants<Char>::digit_0 + unsigned(C)));
                return ret;
            } else {
                constexpr auto prevDigits = Color(unsigned(C) / 10);
                constexpr auto myDigit = Color(unsigned(C) % 10);
                return mytype::encode<myDigit>(mytype::encode<prevDigits>(prev));
            }
        }

        template<Color First, Color... Rest, size_t N>
        static constexpr auto addNextColor(Str<N> prev) {
            auto withSemicolon = mytype::append(prev, CharConstants<Char>::semicolon);
            auto next = mytype::encode<First>(withSemicolon);
            if constexpr (sizeof...(Rest) != 0) {
                return mytype::addNextColor<Rest...>(next);
            } else {
                return mytype::append(next, CharConstants<Char>::letter_m);
            }
        }

        template<Color First, Color... Rest>
        static constexpr auto doMakeColor() {
            constexpr Str<3> prefix{{CharConstants<Char>::esc, CharConstants<Char>::squareBracketOpen, 0}};
            constexpr auto start = mytype::encode<First>(prefix);

            if constexpr (sizeof...(Rest) != 0) {
                return mytype::addNextColor<Rest...>(start);
            } else {
                return mytype::append(start, CharConstants<Char>::letter_m);
            }
        }

    private:
        template<Color First, Color... Rest>
        static constexpr auto storage = mytype::doMakeColor<First, Rest...>();

    public:
        template<Color First, Color... Rest>
        static constexpr auto make() {
            return std::basic_string_view<Char>(mytype::storage<First, Rest...>.wrapped);
        }
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(ColorMaker)

    ARGUM_MOD_EXPORTED
    template<Character Char, Color First, Color... Rest>
    constexpr auto basicMakeColor() {
        return BasicColorMaker<Char>::template make<First, Rest...>();
    }

    ARGUM_MOD_EXPORTED
    template<Color First, Color... Rest>
    constexpr auto makeColor() {
        return BasicColorMaker<char>::template make<First, Rest...>();
    }

    ARGUM_MOD_EXPORTED
    template<Color First, Color... Rest>
    constexpr auto makeWColor() {
        return BasicColorMaker<wchar_t>::template make<First, Rest...>();
    }


    ARGUM_MOD_EXPORTED
    template<Color First, Color... Rest, Character Char>
    auto colorize(std::basic_string_view<Char> str) -> std::basic_string<Char> {
        auto prefix = basicMakeColor<Char, First, Rest...>();
        auto suffix = basicMakeColor<Char, Color::normal>();

        std::basic_string<Char> ret;
        ret.reserve(prefix.size() + str.size() + suffix.size());
        ret = prefix;
        ret += str;
        ret += suffix;
        return ret;
    }


    ARGUM_MOD_EXPORTED
    template<class Char>
    struct BasicColorScheme {
        using StringViewType = std::basic_string_view<Char>;

        StringViewType heading;
        StringViewType progName;
        StringViewType shortOptionInUsage;
        StringViewType longOptionInUsage;
        StringViewType optionArgInUsage;
        StringViewType positionalInUsage;
        StringViewType shortOption;
        StringViewType longOption;
        StringViewType optionArg;
        StringViewType positional;
        StringViewType error;
        StringViewType warning;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(ColorScheme)

    ARGUM_MOD_EXPORTED
    template<Character Char>
    constexpr BasicColorScheme<Char> nullColorScheme {};

    ARGUM_MOD_EXPORTED
    template<Character Char>
    constexpr BasicColorScheme<Char> basicDefaultColorScheme {
            .heading =              basicMakeColor<Char, Color::bold, Color::blue>(),
            .progName =             basicMakeColor<Char, Color::bold, Color::magenta>(),
            .shortOptionInUsage =   basicMakeColor<Char, Color::green>(),
            .longOptionInUsage =    basicMakeColor<Char, Color::cyan>(),
            .optionArgInUsage =     basicMakeColor<Char, Color::yellow>(),
            .positionalInUsage =    basicMakeColor<Char, Color::green>(),
            .shortOption =          basicMakeColor<Char, Color::bold, Color::green>(),
            .longOption =           basicMakeColor<Char, Color::bold, Color::cyan>(),
            .optionArg =            basicMakeColor<Char, Color::bold, Color::yellow>(),
            .positional =           basicMakeColor<Char, Color::bold, Color::green>(),
            .error =                basicMakeColor<Char, Color::bold, Color::red>(),
            .warning =              basicMakeColor<Char, Color::bold, Color::yellow>(),
    };

    ARGUM_MOD_EXPORTED
    constexpr auto defaultColorScheme() -> const ColorScheme & { return basicDefaultColorScheme<char>; }

    ARGUM_MOD_EXPORTED
    constexpr auto defaultWColorScheme() -> const WColorScheme & { return basicDefaultColorScheme<wchar_t>; }

    ARGUM_MOD_EXPORTED
    template<Character Char>
    class BasicColorizer {
    public:
        using CharType = Char;
        using StringViewType = std::basic_string_view<CharType>;
        using StringType = std::basic_string<CharType>;
        using Scheme = BasicColorScheme<Char>;

    public:
        constexpr BasicColorizer() = default;

        constexpr BasicColorizer(const Scheme & scheme):
            m_scheme(&scheme)
        {}

        
        auto heading(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->heading);
        }

        auto progName(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->progName);
        }

        auto shortOptionInUsage(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->shortOptionInUsage);
        }

        auto longOptionInUsage(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->longOptionInUsage);
        }

        auto optionArgInUsage(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->optionArgInUsage);
        }

        auto positionalInUsage(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->positionalInUsage);
        }

        auto shortOption(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->shortOption);
        }

        auto longOption(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->longOption);
        }

        auto optionArg(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->optionArg);
        }

        auto positional(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->positional);
        }

        auto error(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->error);
        }

        auto warning(StringViewType str) const -> StringType {
            return this->colorize(str, m_scheme->warning);
        }
    private:
        static auto colorize(StringViewType str, StringViewType prefix) -> StringType {
            if (!prefix.size())
                return StringType(str);

            auto suffix = basicMakeColor<Char, Color::normal>();

            StringType ret;
            ret.reserve(prefix.size() + str.size() + suffix.size());
            ret = prefix;
            ret += str;
            ret += suffix;
            return ret;
        }

    private:
        const Scheme * m_scheme = &nullColorScheme<Char>;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(Colorizer)

    ARGUM_MOD_EXPORTED 
    inline constexpr auto defaultColorizer() -> Colorizer {
        return {basicDefaultColorScheme<char>};
    }
    ARGUM_MOD_EXPORTED 
    inline constexpr auto defaultWColorizer() -> WColorizer {
        return {basicDefaultColorScheme<wchar_t>};
    }
}

#endif


namespace Argum {

    template<class Char> class BasicOption;
    template<class Char> class BasicPositional;
    template<class Char> class BasicParser;

    ARGUM_MOD_EXPORTED
    template<class Char>
    class BasicHelpFormatter {
    public:
        using CharType = Char;
        using StringViewType = std::basic_string_view<CharType>;
        using StringType = std::basic_string<CharType>;
        using Option = BasicOption<CharType>;
        using Positional = BasicPositional<CharType>;
        using Colorizer = BasicColorizer<CharType>;

        struct Layout {
            unsigned width = std::numeric_limits<unsigned>::max();
            unsigned helpLeadingGap = 2;
            unsigned helpNameMaxWidth = 20;
            unsigned helpDescriptionGap = 2;
        };
        static constexpr Layout defaultLayout = {};

        struct SubCommandMark {
            size_t positionalIdx = size_t(-1);
            size_t optionIdx = size_t(-1);
        };
    private:
        using CharConstants = Argum::CharConstants<CharType>;
        using Messages = Argum::Messages<CharType>;

    public:
        BasicHelpFormatter(const BasicParser<Char> & parser, StringViewType progName, Layout layout = defaultLayout):
            m_progName(progName),
            m_parser(parser),
            m_layout(layout) {
            
            if (m_layout.helpNameMaxWidth == 0)
                m_layout.helpNameMaxWidth = 1;
            if (m_layout.width <= m_layout.helpLeadingGap + m_layout.helpNameMaxWidth + m_layout.helpDescriptionGap)
                m_layout.width = m_layout.helpLeadingGap + m_layout.helpNameMaxWidth + m_layout.helpDescriptionGap + 1;
        }

        auto formatUsage(const Colorizer & colorizer = {}) const -> StringType {
            return this->formatUsage(std::nullopt, colorizer);
        }

        auto formatUsage(const std::optional<StringType> & subCommand,
                         const Colorizer & colorizer = {}) const -> StringType {
            constexpr auto space = CharConstants::space;
            return wordWrap(colorizer.heading(Messages::usageStart()).
                                append(colorizer.progName(this->m_progName)).
                                append({space}).
                                append(this->formatSyntax(subCommand, colorizer)), 
                            m_layout.width, m_layout.helpLeadingGap);
        }

        auto formatHelp(const Colorizer & colorizer = {}) const -> StringType {
            return this->formatHelp(std::nullopt, colorizer);
        }

        auto formatHelp(const std::optional<StringType> & subCommand,
                        const Colorizer & colorizer = {}) const -> StringType {

            if (subCommand && this->m_parser.subCommandMark().positionalIdx == size_t(-1))
                ARGUM_INVALID_ARGUMENT("subcommand must be defined to use this function with non null subcommand");

            constexpr auto endl = CharConstants::endl;

            StringType ret;

            auto helpContent = this->calculateHelpContent(bool(subCommand), colorizer);
            if (helpContent.maxNameLen > m_layout.helpNameMaxWidth)
                helpContent.maxNameLen = m_layout.helpNameMaxWidth;

            if (!helpContent.positionalItems.empty()) {
                ret.append(wordWrap(colorizer.heading(Messages::positionalHeader()), m_layout.width, m_layout.helpLeadingGap));
                for(auto & [name, desc]: helpContent.positionalItems) {
                    ret.append({endl}).append(this->formatItemHelp(name, desc, helpContent.maxNameLen));
                }
                ret.append(2, endl);
            }

            if (!helpContent.optionItems.empty()) {
                ret.append(wordWrap(colorizer.heading(Messages::optionsHeader()), m_layout.width, m_layout.helpLeadingGap));
                for(auto & [name, desc]: helpContent.optionItems) {
                    ret.append({endl}).append(this->formatItemHelp(name, desc, helpContent.maxNameLen));
                }
                ret.append(2, endl);
            }
            return ret;
        }

        auto formatSyntax(const Colorizer & colorizer = {}) const -> StringType {

            return this->formatSyntax(std::nullopt, colorizer);
        }

        auto formatSyntax(const std::optional<StringType> & subCommand,
                          const Colorizer & colorizer = {}) const -> StringType {

            auto subCommandMark = this->m_parser.subCommandMark();
            if (subCommand && subCommandMark.positionalIdx == size_t(-1))
                ARGUM_INVALID_ARGUMENT("subcommand must be added to use this function with non null subcommand");

            constexpr auto space = CharConstants::space;

            auto getSyntax = [&](auto & obj) {
                return obj.formatSyntax(this->m_parser, colorizer);
            };

            StringType ret = this->appendSyntax(
                    join(this->optionsBegin(false), this->optionsEnd(false), space, getSyntax),
                    join(this->positionalsBegin(false), this->positionalsEnd(false), space, getSyntax)
            );
            if (subCommand) {
                ret = this->appendSyntax(std::move(ret), *subCommand); 
                ret = this->appendSyntax(std::move(ret), join(this->optionsBegin(true), this->optionsEnd(true), space, getSyntax));
                ret = this->appendSyntax(std::move(ret), join(this->positionalsBegin(true), this->positionalsEnd(true), space, getSyntax));
            } else if (subCommandMark.positionalIdx != size_t(-1)) {
                ret = this->appendSyntax(std::move(ret), getSyntax(this->m_parser.positionals()[subCommandMark.positionalIdx]));
            }

            return ret;
        }
        

        struct HelpContent {
            unsigned maxNameLen = 0;
            std::vector<std::pair<StringType, StringType>> optionItems;
            std::vector<std::pair<StringType, StringType>> positionalItems;
        };
        auto calculateHelpContent(bool forSubCommand, 
                                  const Colorizer & colorizer = {}) const -> HelpContent {
            
            HelpContent ret;
            auto subCommandMark = this->m_parser.subCommandMark();

            size_t positionalsSize;
            if (forSubCommand || subCommandMark.positionalIdx == size_t(-1))
                positionalsSize = this->m_parser.positionals().size();
            else 
                positionalsSize = subCommandMark.positionalIdx + 1;
            for(size_t i = 0; i < positionalsSize; ++ i) {
                if (forSubCommand && i == subCommandMark.positionalIdx)
                    continue;
                auto & pos = this->m_parser.positionals()[i];
                auto name = pos.formatHelpName(this->m_parser, colorizer);
                auto length = stringWidth(name);
                if (length > ret.maxNameLen)
                    ret.maxNameLen = length;
                ret.positionalItems.emplace_back(std::move(name), pos.formatHelpDescription());
            }
            std::for_each(this->m_parser.options().begin(), this->optionsEnd(forSubCommand), [&](auto & opt){
                auto name = opt.formatHelpName(this->m_parser, colorizer);
                auto length = stringWidth(name);
                if (length > ret.maxNameLen)
                    ret.maxNameLen = length;
                ret.optionItems.emplace_back(std::move(name), opt.formatHelpDescription());
            });
            return ret;
        }

        auto formatItemHelp(StringViewType name, 
                            StringViewType description,
                            unsigned maxNameLen) const -> StringType {
            constexpr auto space = CharConstants::space;
            constexpr auto endl = CharConstants::endl;

            auto descColumnOffset = this->m_layout.helpLeadingGap + maxNameLen + this->m_layout.helpDescriptionGap;

            StringType ret = wordWrap(StringType(this->m_layout.helpLeadingGap, space).append(name), this->m_layout.width, this->m_layout.helpLeadingGap);
            auto lastEndlPos = ret.rfind(endl);
            auto lastLineLen = stringWidth(StringViewType(ret.c_str() + (lastEndlPos + 1), ret.size() - (lastEndlPos + 1)));

            if (lastLineLen > maxNameLen + this->m_layout.helpLeadingGap) {
                ret += endl;
                ret.append(descColumnOffset, space);
            } else {
                ret.append(descColumnOffset - lastLineLen, space);
            }

            ret.append(wordWrap(description, this->m_layout.width, descColumnOffset, descColumnOffset));

            return ret;
        }

        auto optionsBegin(bool forSubCommand) const {
            if (!forSubCommand)
                return  this->m_parser.options().begin();
            auto subCommandMark = this->m_parser.subCommandMark();
            if (subCommandMark.optionIdx != size_t(-1))
                return this->m_parser.options().begin() + subCommandMark.optionIdx;
            return this->m_parser.options().end();
        }
        auto optionsEnd(bool forSubCommand) const {
            auto subCommandMark = this->m_parser.subCommandMark();
            if (forSubCommand || subCommandMark.optionIdx == size_t(-1))
                return this->m_parser.options().end();
            return this->m_parser.options().begin() + subCommandMark.optionIdx;
        }
        auto positionalsBegin(bool forSubCommand) const {
            if (!forSubCommand)
                return this->m_parser.positionals().begin();
            auto subCommandMark = this->m_parser.subCommandMark();
            if (subCommandMark.positionalIdx != size_t(-1))
                return this->m_parser.positionals().begin() + subCommandMark.positionalIdx + 1;
            return this->m_parser.positionals().end();
        }
        auto positionalsEnd(bool forSubCommand) const {
            auto subCommandMark = this->m_parser.subCommandMark();
            if (forSubCommand || subCommandMark.positionalIdx == size_t(-1))
                return this->m_parser.positionals().end();
            return this->m_parser.positionals().begin() + subCommandMark.positionalIdx;
        }

        static auto appendSyntax(StringType base, StringType addend) -> StringType {
            
            constexpr auto space = CharConstants::space;

            if (!base.empty()) {
                if (!addend.empty()) {
                    base += space;
                    base += std::move(addend);
                    return base;
                }
                return base;
            } 
            
            return addend;
        }
    private:
        StringType m_progName;
        const BasicParser<Char> & m_parser;
        Layout m_layout;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(HelpFormatter)

}


#endif




namespace Argum {

    template<class Char> class BasicParser;

    template<class Char>
    class BasicOption {
        friend class BasicParser<Char>;

    public:
        using CharType = Char;
        using StringType = std::basic_string<Char>;
        using StringViewType = std::basic_string_view<Char>;
        using OptionNames = BasicOptionNames<Char>;
        using ArgumentKind = OptionArgumentKind;

    private:
        using CharConstants = Argum::CharConstants<CharType>;
        using Messages = Argum::Messages<CharType>;
        using Colorizer = BasicColorizer<CharType>;
        using HandlerReturnType = ARGUM_EXPECTED(CharType, void);

        using NoArgHandler = std::function<HandlerReturnType ()>; 
        using OptArgHandler = std::function<HandlerReturnType (std::optional<StringViewType>)>; 
        using ReqArgHandler = std::function<HandlerReturnType (StringViewType)>; 

        template<class HandlerType>
        static constexpr auto argumentKindOf() -> ArgumentKind {
            if constexpr (std::is_same_v<HandlerType, NoArgHandler>)
                return ArgumentKind::None;
            else if constexpr (std::is_same_v<HandlerType, OptArgHandler>) 
                return ArgumentKind::Optional;
            else if constexpr (std::is_same_v<HandlerType, ReqArgHandler>) 
                return ArgumentKind::Required;
        }

        auto canHaveArgument() const -> bool {
            return !std::holds_alternative<NoArgHandler>(m_handler);
        }

    public:
        using Handler = std::variant<
            NoArgHandler,
            OptArgHandler,
            ReqArgHandler
        >;
        
        BasicOption(OptionNames names) :
            m_names(std::move(names)) {
        }

        template<class... Args>
        requires(std::is_constructible_v<OptionNames, Args &&...>)
        BasicOption(Args && ...args) :
            m_names(std::forward<Args>(args)...) {
        }

        template<class H>
        auto handler(H && h) & -> BasicOption &
        requires(IsHandler<decltype(h), CharType, void> || 
                 IsHandler<decltype(h), CharType, void, std::optional<StringViewType>> ||
                 IsHandler<decltype(h), CharType, void, StringViewType>) {

            using InHandlerType = std::remove_cvref_t<decltype(h)>;

            if constexpr (std::is_invocable_v<InHandlerType>) {
                this->m_handler.template emplace<NoArgHandler>(adaptHandler<CharType, H, void>(std::forward<H>(h)));
            } else if constexpr (std::is_invocable_v<InHandlerType, std::optional<StringViewType>>)
                this->m_handler.template emplace<OptArgHandler>(adaptHandler<CharType, H, void, std::optional<StringViewType>>(std::forward<H>(h)));
            else {
                this->m_handler.template emplace<ReqArgHandler>(adaptHandler<CharType, H, void, StringViewType>(std::forward<H>(h)));
            }
            return *this;
        }

        template<class H>
        auto handler(H && h) && -> BasicOption &&
        requires(IsHandler<decltype(h), CharType, void> || 
                 IsHandler<decltype(h), CharType, void, std::optional<StringViewType>> ||
                 IsHandler<decltype(h), CharType, void, StringViewType>) {
            return std::move(static_cast<BasicOption *>(this)->handler(std::forward<H>(h)));
        }

        auto occurs(Quantifier r) & -> BasicOption & {
            this->m_occurs = r;
            return *this;
        }
        auto occurs(Quantifier r) && -> BasicOption && {
            return std::move(static_cast<BasicOption *>(this)->occurs(r));
        }

        auto argName(StringViewType n) & -> BasicOption & {
            this->m_argName = n;
            return *this;
        }
        auto argName(StringViewType n) && -> BasicOption && {
            return std::move(static_cast<BasicOption *>(this)->argName(n));
        }

        auto requireAttachedArgument(bool val) & -> BasicOption & {
            this->m_requireAttachedArgument = val;
            return *this;
        }
        auto requireAttachedArgument(bool val) && -> BasicOption && {
            return std::move(static_cast<BasicOption *>(this)->requireAttachedArgument(val));
        }

        auto help(StringViewType str) & -> BasicOption & {
            this->m_description = str;
            return *this;
        }
        auto help(StringViewType str) && -> BasicOption && {
            return std::move(static_cast<BasicOption *>(this)->help(str));
        }

        auto formatSyntax(const BasicParser<CharType> & parser,
                          const Colorizer & colorizer = {}) const -> StringType {
            constexpr auto space = CharConstants::space;
            constexpr auto brop = CharConstants::squareBracketOpen;
            constexpr auto brcl = CharConstants::squareBracketClose;
            
            StringType ret;

            if (this->m_occurs.min() == 0)
                ret += brop;
            auto & mainName = this->m_names.main();
            bool isLong = parser.isOptionNameLong(mainName);
            StringType nameAndArg = isLong ? colorizer.longOptionInUsage(mainName) : colorizer.shortOptionInUsage(mainName);
            nameAndArg.append(this->formatArgSyntax(isLong, true, colorizer));
            ret.append(nameAndArg);
            unsigned idx = 1;
            for (; idx < this->m_occurs.min(); ++idx) {
                ret.append({space}).append(nameAndArg);
            }
            if (this->m_occurs.max() != Quantifier::infinity && 
                idx < this->m_occurs.max()) {
                
                ret.append({space, brop}).append(nameAndArg);
                for (++idx; idx < this->m_occurs.max(); ++idx)
                    ret.append({space}).append(nameAndArg);
                ret += brcl;
            }
            if (this->m_occurs.min() == 0)
                ret += brcl;

            return ret;
        }

        auto formatArgSyntax(bool forLongName, bool forUsage, const Colorizer & colorizer = {}) const -> StringType {
            constexpr auto space = CharConstants::space;
            constexpr auto brop = CharConstants::squareBracketOpen;
            constexpr auto brcl = CharConstants::squareBracketClose;
            constexpr auto eq = CharConstants::assignment;

            StringType ret;
            std::visit([&](const auto & handler) {
                using HandlerType = std::remove_cvref_t<decltype(handler)>;
                constexpr auto argumentKind = BasicOption::template argumentKindOf<HandlerType>();
                auto colorArg = forUsage ? colorizer.optionArgInUsage(this->m_argName) : colorizer.optionArg(this->m_argName);
                if constexpr (argumentKind == ArgumentKind::Optional) {
                    if (this->m_requireAttachedArgument)
                        if (forLongName)
                            ret.append({brop, eq});
                        else
                            ret.append({brop});
                    else
                        ret.append({space, brop});
                    ret.append(colorArg).append({brcl});
                } else if constexpr (argumentKind == ArgumentKind::Required)  {
                    if (this->m_requireAttachedArgument) {
                        if (forLongName)
                            ret.append({eq});
                    } else {
                        ret.append({space});
                    }
                    ret.append(colorArg);
                }
            }, this->m_handler);
            return ret;
        }

        auto formatHelpName(const BasicParser<CharType> & parser, 
                            const Colorizer & colorizer = {}) const -> StringType {

            auto handleName = [&](const StringType & name) {
                auto isLong = parser.isOptionNameLong(name);
                StringType ret = isLong ? colorizer.longOption(name) : colorizer.shortOption(name);
                ret += this->formatArgSyntax(isLong, false, colorizer);
                return ret;
            };

            auto & all = this->m_names.all();
            auto first = all.begin();
            auto last = all.end();
        
            StringType ret = handleName(*first);
            ++first;
            std::for_each(first, last, [&](const auto & name) {
                ret.append(Messages::listJoiner()).append(handleName(name));
            });

            return ret;
        }

        auto formatHelpDescription() const -> const StringType & {
            return this->m_description;
        }
    private:
        OptionNames m_names;
        Handler m_handler = []() -> ARGUM_EXPECTED(CharType, void) { return ARGUM_VOID_SUCCESS; };
        Quantifier m_occurs = zeroOrMoreTimes;

        StringType m_argName = Messages::defaultArgName();
        StringType m_description;
        bool m_requireAttachedArgument = false;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(Option)

    template<class Char>
    class BasicPositional {
        friend class BasicParser<Char>;
    
    public:
        using CharType = Char;
        using StringType = std::basic_string<Char>;
        using StringViewType = std::basic_string_view<Char>;
        using Handler = std::function<ARGUM_EXPECTED(CharType, void) (StringViewType)>;

    private:
        using CharConstants = Argum::CharConstants<CharType>;
        using Colorizer = BasicColorizer<CharType>;
        using HandlerReturnType = ARGUM_EXPECTED(CharType, void);

    public:
        BasicPositional(StringViewType name) :
            m_name(std::move(name)) {
        }

        template<class H>
        auto handler(H && h) & -> BasicPositional & 
        requires(IsHandler<decltype(h), CharType, void, StringViewType>) {
            this->m_handler = adaptHandler<CharType, H, void, StringViewType>(std::forward<H>(h));
            return *this;
        }
        template<class H>
        auto handler(H && h) && -> BasicPositional && 
        requires(IsHandler<decltype(h), CharType, void, StringViewType>) {
            return std::move(static_cast<BasicPositional *>(this)->handler(std::forward<H>(h)));
        }

        auto occurs(Quantifier r) & -> BasicPositional & {
            this->m_occurs = r;
            return *this;
        }
        auto occurs(Quantifier r) && -> BasicPositional && {
            return std::move(static_cast<BasicPositional *>(this)->occurs(r));
        }

        auto help(StringViewType str) & -> BasicPositional &{
            this->m_description = str;
            return *this;
        }
        auto help(StringViewType str) && -> BasicPositional && {
            return std::move(static_cast<BasicPositional *>(this)->help(str));
        }

        auto formatSyntax(const BasicParser<CharType> & /*parser*/,
                          const Colorizer & colorizer = {}) const -> StringType {
            constexpr auto space = CharConstants::space;
            constexpr auto brop = CharConstants::squareBracketOpen;
            constexpr auto brcl = CharConstants::squareBracketClose;
            constexpr auto ellipsis = CharConstants::ellipsis;

            StringType ret;

            StringType colorizedName = colorizer.positionalInUsage(this->m_name);

            if (this->m_occurs.min() == 0)
                ret += brop;
            ret += colorizedName;
            unsigned idx = 1;
            for (; idx < this->m_occurs.min(); ++idx) {
                ret.append({space}).append(colorizedName);
            }
            if (idx < this->m_occurs.max()) {
                ret.append({space, brop}).append(colorizedName);
                if (this->m_occurs.max() != Quantifier::infinity) {
                    for (++idx; idx < this->m_occurs.max(); ++idx)
                        ret.append({space}).append(colorizedName);
                } else {
                    ret.append({space}).append(colorizer.positionalInUsage(ellipsis));
                }
                ret += brcl;
            }
            if (this->m_occurs.min() == 0)
                ret += brcl;

            return ret;
        }

        auto formatHelpName(const BasicParser<CharType> & /*parser*/,
                          const Colorizer & colorizer = {}) const -> StringType {
            return colorizer.positional(this->m_name);
        }

        auto formatHelpDescription() const -> const StringType & {
            return this->m_description;
        }
    private:
        StringType m_name;
        Handler m_handler = [](StringViewType) -> ARGUM_EXPECTED(CharType, void) { return ARGUM_VOID_SUCCESS; };
        Quantifier m_occurs = once;
        StringType m_description;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(Positional)

    template<class Char>
    class BasicParser {

    public:
        using CharType = Char;
        using StringType = std::basic_string<Char>;
        using StringViewType = std::basic_string_view<Char>;
        using OptionNames = BasicOptionNames<Char>;
        using ParsingException = BasicParsingException<Char>;
        using Option = BasicOption<Char>;
        using Positional = BasicPositional<Char>;
        using Settings = typename BasicTokenizer<Char>::Settings;
        using HelpFormatter = BasicHelpFormatter<CharType>;
        using SubCommandMark = typename HelpFormatter::SubCommandMark;
        using Colorizer = BasicColorizer<Char>;

    private:
        using CharConstants = Argum::CharConstants<CharType>;
        using Messages = Argum::Messages<CharType>;
        using Tokenizer = BasicTokenizer<CharType>;
        using ValidationData = BasicValidationData<CharType>;

        using ValidatorFunction = std::function<bool (const ValidationData &)>;

    public:
        struct UnrecognizedOption : public ParsingException {
            ARGUM_IMPLEMENT_EXCEPTION(UnrecognizedOption, ParsingException, Error::UnrecognizedOption)

            UnrecognizedOption(StringViewType option_): 
                ParsingException(ErrorCode, format(Messages::unrecognizedOptionError(), option_)),
                option(option_) {
            }
            StringType option;
        };

        struct AmbiguousOption : public ParsingException {
            ARGUM_IMPLEMENT_EXCEPTION(AmbiguousOption, ParsingException, Error::AmbiguousOption)

            AmbiguousOption(StringViewType option_, std::vector<StringType> possibilities_): 
                ParsingException(ErrorCode, 
                                 format(Messages::ambiguousOptionError(), option_, 
                                        join(possibilities_.begin(), possibilities_.end(), Messages::listJoiner()))),
                option(option_),
                possibilities(std::move(possibilities_)) {
            }
            StringType option;
            std::vector<StringType> possibilities;
        };

        struct MissingOptionArgument : public ParsingException {
            ARGUM_IMPLEMENT_EXCEPTION(MissingOptionArgument, ParsingException, Error::MissingOptionArgument)

            MissingOptionArgument(StringViewType option_): 
                ParsingException(ErrorCode, format(Messages::missingOptionArgumentError(), option_)),
                option(option_) {
            }
            StringType option;
        };

        struct ExtraOptionArgument : public ParsingException {
            ARGUM_IMPLEMENT_EXCEPTION(ExtraOptionArgument, ParsingException, Error::ExtraOptionArgument)

            ExtraOptionArgument(StringViewType option_): 
                ParsingException(ErrorCode, format(Messages::extraOptionArgumentError(), option_)),
                option(option_) {
            }
            StringType option;
        };

        struct ExtraPositional : public ParsingException {
            ARGUM_IMPLEMENT_EXCEPTION(ExtraPositional, ParsingException, Error::ExtraPositional)

            ExtraPositional(StringViewType value_): 
                ParsingException(ErrorCode, format(Messages::extraPositionalError(), value_)),
                value(value_) {
            }
            StringType value;
        };
        
        struct ValidationError : public ParsingException {
            ARGUM_IMPLEMENT_EXCEPTION(ValidationError, ParsingException, Error::ValidationError)

            ValidationError(StringViewType message):
                ParsingException(ErrorCode, format(Messages::validationError(), message)) {
            }
            template<DescribableParserValidator<CharType> Validator>
            ValidationError(Validator validator):
                ParsingException(ErrorCode, format(Messages::validationError(), describe(validator))) {
            }
        };

    public:
        BasicParser() = default;

        BasicParser(Settings settings): m_tokenizer(settings) {
        }

        auto add(Option option) -> void {

            auto & added = this->m_options.emplace_back(std::move(option));
            this->m_tokenizer.add(this->m_options.back().m_names);
            if (added.m_occurs.min() > 0)
                addValidator(optionOccursAtLeast(added.m_names.main(), added.m_occurs.min()));
            ++m_updateCount;
        }
        
        auto add(Positional positional) {
            //for expected number of positionals this is faster than maintaining a map
            for(auto & existing: this->m_positionals) {
                if (existing.m_name == positional.m_name)
                    ARGUM_INVALID_ARGUMENT("duplicate positional name");
            }
            this->m_positionals.emplace_back(std::move(positional));
            ++m_updateCount;
        }

        auto addSubCommand(Positional positional) {
            this->add(std::move(positional));
            this->m_subCommandMark = {this->m_positionals.size() - 1,  this->m_options.size()};
        }

        template<ParserValidator<CharType> Validator>
        auto addValidator(Validator v, StringViewType description) -> void {
            m_validators.emplace_back(std::move(v), description);
        }

        template<DescribableParserValidator<CharType> Validator>
        auto addValidator(Validator v) -> void  {
            auto desc = describe(v);
            m_validators.emplace_back(std::move(v), std::move(desc));
        }

        auto parse(int argc, CharType ** argv) const -> ARGUM_EXPECTED(CharType, void) {
            return this->parse(makeArgSpan<CharType>(argc, argv));
        }

        template<ArgRange<CharType> Args>
        auto parse(const Args & args) const -> ARGUM_EXPECTED(CharType, void) {
            
            return this->parse(std::begin(args), std::end(args));
        }

        template<ArgIterator<CharType> It>
        auto parse(It argFirst, It argLast) const -> ARGUM_EXPECTED(CharType, void) {
            
            using ReturnType = ARGUM_EXPECTED(CharType, void);
            ParsingState parsingState(*this);

            return ReturnType(parsingState.parse(argFirst, argLast, /*stopOnUnknown=*/false));
        }

        auto parseUntilUnknown(int argc, CharType ** argv) const -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            return this->parseUntilUnknown(makeArgSpan<CharType>(argc, argv));
        }

        template<ArgRange<CharType> Args>
        auto parseUntilUnknown(const Args & args) const -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            
            return this->parseUntilUnknown(std::begin(args), std::end(args));
        }

        template<ArgIterator<CharType> It>
        auto parseUntilUnknown(It argFirst, It argLast) const -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            
            ParsingState parsingState(*this);

            return parsingState.parse(argFirst, argLast, /*stopOnUnknown=*/true);
        }

        auto options() const -> const std::vector<Option> & {
            return this->m_options;
        }

        auto positionals() const -> const std::vector<Positional> & {
            return this->m_positionals;
        }

        auto subCommandMark() const -> SubCommandMark {
            return this->m_subCommandMark;
        }

        auto formatUsage(StringViewType progName) const -> StringType {
            return this->formatUsage(progName, {}, HelpFormatter::defaultLayout.width, {});
        }

        auto formatUsage(StringViewType progName, unsigned width) const -> StringType {
            return this->formatUsage(progName, {}, width, {});
        }

        auto formatUsage(StringViewType progName, const Colorizer & colorizer) const -> StringType {
            return this->formatUsage(progName, {}, HelpFormatter::defaultLayout.width, colorizer);
        }

        auto formatUsage(StringViewType progName, 
                         unsigned width,
                         const Colorizer & colorizer) const -> StringType {
            return this->formatUsage(progName, {}, width, colorizer);
        }

        auto formatUsage(StringViewType progName, std::optional<StringType> subCommand) const -> StringType {
            return this->formatUsage(progName, subCommand, HelpFormatter::defaultLayout.width, {});
        }

        auto formatUsage(StringViewType progName, std::optional<StringType> subCommand,
                         unsigned width) const -> StringType {
            return this->formatUsage(progName, subCommand, width, {});
        }

        auto formatUsage(StringViewType progName, std::optional<StringType> subCommand,
                         const Colorizer & colorizer) const -> StringType {
            return this->formatUsage(progName, subCommand, HelpFormatter::defaultLayout.width, colorizer);
        }

        auto formatUsage(StringViewType progName, std::optional<StringType> subCommand,
                         unsigned width,
                         const Colorizer & colorizer) const -> StringType {
            typename HelpFormatter::Layout layout;
            layout.width = width;
            return HelpFormatter(*this, progName, layout).formatUsage(subCommand, colorizer);
        }

        auto formatHelp(StringViewType progName) {
            return this->formatHelp(progName, {}, HelpFormatter::defaultLayout.width, {});
        }

        auto formatHelp(StringViewType progName, unsigned width) {
            return this->formatHelp(progName, {}, width, {});
        }

        auto formatHelp(StringViewType progName, const Colorizer & colorizer) {
            return this->formatHelp(progName, {}, HelpFormatter::defaultLayout.width, colorizer);
        }

        auto formatHelp(StringViewType progName, unsigned width, const Colorizer & colorizer) {
            return this->formatHelp(progName, {}, width, colorizer);
        }

        auto formatHelp(StringViewType progName, std::optional<StringType> subCommand) {
            return this->formatHelp(progName, subCommand, HelpFormatter::defaultLayout.width, {});
        }

        auto formatHelp(StringViewType progName, std::optional<StringType> subCommand,
                        unsigned width) {
            return this->formatHelp(progName, subCommand, width, {});
        }

        auto formatHelp(StringViewType progName, std::optional<StringType> subCommand,
                        const Colorizer & colorizer) {
            return this->formatHelp(progName, subCommand, HelpFormatter::defaultLayout.width, colorizer);
        }

        auto formatHelp(StringViewType progName, std::optional<StringType> subCommand,
                        unsigned width,
                        const Colorizer & colorizer) const -> StringType {
            typename HelpFormatter::Layout layout;
            layout.width = width;           
            HelpFormatter formatter(*this, progName, layout);
            StringType ret = formatter.formatUsage(subCommand, colorizer);
            ret.append(2, CharConstants::endl);
            ret.append(formatter.formatHelp(subCommand, colorizer));
            return ret;
        }

        auto isOptionNameLong(const StringViewType name) const -> bool {
            return this->m_tokenizer.isOptionNameLong(name);
        }

    private:
        class ParsingState {
        public:
            ParsingState(const BasicParser & owner): 
                m_owner(owner),
                m_updateCountAtLastRecalc(owner.m_updateCount - 1) {
            }

            template<ArgIterator<CharType> It>
            auto parse(It argFirst, It argLast, bool stopOnUnknown) -> ARGUM_EXPECTED(CharType, std::vector<StringType>) {
            
                ARGUM_CHECK_RESULT(auto ret, m_owner.m_tokenizer.tokenize(argFirst, argLast, [&](auto && token) -> ARGUM_EXPECTED(CharType, typename Tokenizer::TokenResult) {

                    using TokenType = std::remove_cvref_t<decltype(token)>;

                    if constexpr (std::is_same_v<TokenType, typename Tokenizer::OptionToken>) {

                        ARGUM_PROPAGATE_ERROR(resetOption(token.idx, token.usedName, token.argument));
                        return Tokenizer::Continue;

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::OptionStopToken>) {

                        ARGUM_PROPAGATE_ERROR(completeOption());
                        return Tokenizer::Continue;

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::ArgumentToken>) {

                        ARGUM_CHECK_RESULT(auto result, handlePositional(token.value, argFirst + token.argIdx, argLast));
                        if (!result) {
                            if (stopOnUnknown)
                                return Tokenizer::StopBefore;
                            ARGUM_THROW(ExtraPositional, token.value);
                        }
                        return Tokenizer::Continue;

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::UnknownOptionToken>) {

                        ARGUM_PROPAGATE_ERROR(completeOption());
                        if (stopOnUnknown)
                            return Tokenizer::StopBefore;
                        ARGUM_THROW(UnrecognizedOption, token.name);

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::AmbiguousOptionToken>) {

                        ARGUM_PROPAGATE_ERROR(completeOption());
                        ARGUM_THROW(AmbiguousOption, token.name, std::move(token.possibilities));
                    } 
                }));
                ARGUM_PROPAGATE_ERROR(completeOption());
                ARGUM_PROPAGATE_ERROR(validate());
                return ret;
            }

        private:
            auto resetOption(unsigned index, StringViewType name, const std::optional<StringViewType> & argument) -> ARGUM_EXPECTED(CharType, void) {
                ARGUM_PROPAGATE_ERROR(completeOption());
                m_optionName = std::move(name);
                m_optionArgument = argument;
                m_optionIndex = int(index);
                return ARGUM_VOID_SUCCESS;
            }

            auto completeOption() -> ARGUM_EXPECTED(CharType, void) {
                if (m_optionIndex < 0)
                    return ARGUM_VOID_SUCCESS;

                auto & option = m_owner.m_options[unsigned(m_optionIndex)];
                ARGUM_PROPAGATE_ERROR(validateOptionMax(option));
                ARGUM_PROPAGATE_ERROR(std::visit([&](const auto & handler) -> ARGUM_EXPECTED(CharType, void) {
                    using HandlerType = std::remove_cvref_t<decltype(handler)>;
                    constexpr auto argumentKind = Option::template argumentKindOf<HandlerType>();
                    if constexpr (argumentKind == OptionArgumentKind::None) {
                        if (m_optionArgument)
                            ARGUM_THROW(ExtraOptionArgument, m_optionName);
                        ARGUM_PROPAGATE_ERROR(handler());
                    } else if constexpr (argumentKind == OptionArgumentKind::Optional) {
                        ARGUM_PROPAGATE_ERROR(handler(m_optionArgument));
                    } else {
                        if (!m_optionArgument)
                            ARGUM_THROW(MissingOptionArgument, m_optionName);
                        ARGUM_PROPAGATE_ERROR(handler(*m_optionArgument));
                    }
                    return ARGUM_VOID_SUCCESS;
                }, option.m_handler));
                m_optionIndex = -1;
                m_optionArgument.reset();
                return ARGUM_VOID_SUCCESS;
            }

            auto completeOptionUsingArgument(StringViewType argument) -> ARGUM_EXPECTED(CharType, bool) {

                if (m_optionIndex < 0)
                    return false;

                auto & option = m_owner.m_options[unsigned(m_optionIndex)];
                ARGUM_PROPAGATE_ERROR(validateOptionMax(option));
                bool requireAttachedArgument = option.m_requireAttachedArgument;
                ARGUM_CHECK_RESULT(auto ret, std::visit([&](const auto & handler) -> ARGUM_EXPECTED(CharType, bool) {
                        using HandlerType = std::remove_cvref_t<decltype(handler)>;
                        constexpr auto argumentKind = Option::template argumentKindOf<HandlerType>();
                        if constexpr (argumentKind == OptionArgumentKind::None) {
                            if (m_optionArgument)
                                ARGUM_THROW(ExtraOptionArgument, m_optionName);
                            ARGUM_PROPAGATE_ERROR(handler());
                            return false;
                        } else if constexpr (argumentKind == OptionArgumentKind::Optional) {
                            if (requireAttachedArgument || m_optionArgument) {
                                ARGUM_PROPAGATE_ERROR(handler(m_optionArgument));
                                return false;
                            } else {
                                ARGUM_PROPAGATE_ERROR(handler(argument));
                                return true;
                            }
                        } else {
                            if (m_optionArgument) {
                                ARGUM_PROPAGATE_ERROR(handler(*m_optionArgument));
                                return false;
                            } else if (requireAttachedArgument) {
                                ARGUM_THROW(MissingOptionArgument, m_optionName);
                            } else {
                                ARGUM_PROPAGATE_ERROR(handler(argument));
                                return true;
                            }
                        }
                    }, option.m_handler));
                m_optionIndex = -1;
                m_optionArgument.reset();
                return ret;
            }

            auto validateOptionMax(const Option & option) -> ARGUM_EXPECTED(CharType, void) {
                auto & name = option.m_names.main();
                ++m_validationData.optionCount(name);
                auto validator = optionOccursAtMost(name, option.m_occurs.max());
                if (!validator(m_validationData)) {
                    ARGUM_THROW(ValidationError, validator);
                }
                return ARGUM_VOID_SUCCESS;
            }

            template<ArgIterator<CharType> It>
            auto handlePositional(StringViewType value, It remainingArgFirst, It argLast) -> ARGUM_EXPECTED(CharType, bool) {
                ARGUM_CHECK_RESULT(auto completed, completeOptionUsingArgument(value));
                if (completed)
                    return true;

                calculateRemainingPositionals(remainingArgFirst, argLast);
                
                const Positional * positional = nullptr;
                if (m_positionalIndex >= 0) {
                    if (unsigned(m_positionalIndex) >= m_positionalSizes.size())
                        return false;

                    auto & current = m_owner.m_positionals[unsigned(m_positionalIndex)];
                    if (m_positionalSizes[unsigned(m_positionalIndex)] > m_validationData.positionalCount(current.m_name))
                        positional = &current;

                }

                if (!positional) {
                    auto next = std::find_if(m_positionalSizes.begin() + (m_positionalIndex + 1), m_positionalSizes.end(), [](const unsigned val) {
                        return val > 0;
                    });
                    m_positionalIndex = int(next - m_positionalSizes.begin());
                    if (unsigned(m_positionalIndex) >= m_positionalSizes.size())
                        return false;
                    positional = &m_owner.m_positionals[unsigned(m_positionalIndex)];
                }
                
                auto & count = m_validationData.positionalCount(positional->m_name);
                ARGUM_PROPAGATE_ERROR(positional->m_handler(value));
                ++count;
                return true;
            }

            template<ArgIterator<CharType> It>
            auto calculateRemainingPositionals(It remainingArgFirst, It argLast)  {

                if (m_updateCountAtLastRecalc == m_owner.m_updateCount)
                    return;
                
                //1. Count remaining positional arguments
                unsigned remainingPositionalCount = countRemainingPositionals(remainingArgFirst, argLast);
                
                //2. Build the partitioner for positional ranges
                Partitioner<unsigned> partitioner;
                
                auto fillStartIndex = m_positionalIndex + 1;
                if (m_positionalIndex >= 0) {
                    auto & positional = m_owner.m_positionals[unsigned(m_positionalIndex)];
                    auto count = m_validationData.positionalCount(positional.m_name);
                    if (positional.m_occurs.max() > count) {
                        remainingPositionalCount += count; //account for already processed
                        partitioner.addRange(positional.m_occurs.min(), positional.m_occurs.max());
                        --fillStartIndex;
                    }
                }
                std::for_each(m_owner.m_positionals.begin() + (m_positionalIndex + 1), m_owner.m_positionals.end(),
                                [&] (const Positional & positional) {
                    partitioner.addRange(positional.m_occurs.min(), positional.m_occurs.max());
                });

                //3. Partition the range
                unsigned maxRemainingPositionalCount = std::max(remainingPositionalCount, partitioner.minimumSequenceSize());
                auto res = partitioner.partition(maxRemainingPositionalCount);
                ARGUM_ALWAYS_ASSERT(res); //this must be true by construction

                //4. Fill in expected sizes based on regex matches
                m_positionalSizes.resize(m_owner.m_positionals.size());
                std::copy(res->begin(), res->end() - 1, m_positionalSizes.begin() + fillStartIndex);

                this->m_updateCountAtLastRecalc = m_owner.m_updateCount;
            }
            
            template<ArgIterator<CharType> It>
            auto countRemainingPositionals(It remainingArgFirst, It argLast) -> unsigned {
                
                unsigned remainingPositionalCount = 0;
                bool currentOptionExpectsArgument = false;
                (void)m_owner.m_tokenizer.tokenize(remainingArgFirst, argLast, [&](const auto & token) noexcept {
                    using TokenType = std::remove_cvref_t<decltype(token)>;

                    if constexpr (std::is_same_v<TokenType, typename Tokenizer::OptionToken>) {

                        auto & option = m_owner.m_options[token.idx];
                        if (option.canHaveArgument()) {
                            currentOptionExpectsArgument = !token.argument;
                        } else {
                            currentOptionExpectsArgument = false;
                        }

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::OptionStopToken>) {

                        currentOptionExpectsArgument = false;

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::ArgumentToken>) {

                        if (!currentOptionExpectsArgument)
                            ++remainingPositionalCount;
                        else
                            currentOptionExpectsArgument = false;

                    } else if constexpr (std::is_same_v<TokenType, typename Tokenizer::UnknownOptionToken>) {

                        currentOptionExpectsArgument = false;
                    }
                    return Tokenizer::Continue;
                });
                
                return remainingPositionalCount;
            }

            auto validate() -> ARGUM_EXPECTED(CharType, void) {

                //We could use normal validators for this but it is faster to do it manually
                for(auto idx = (m_positionalIndex >= 0 ? unsigned(m_positionalIndex) : 0u); 
                    idx != unsigned(m_owner.m_positionals.size());
                    ++idx) {
                    
                    auto & positional = m_owner.m_positionals[unsigned(idx)];
                    auto validator = positionalOccursAtLeast(positional.m_name, positional.m_occurs.min());
                    if (!validator(m_validationData))
                        ARGUM_THROW(ValidationError, validator);
                }
                
                for(auto & [validator, desc]: m_owner.m_validators) {
                    if (!validator(m_validationData))
                        ARGUM_THROW(ValidationError, desc);
                }

                return ARGUM_VOID_SUCCESS;
            }

        private:
            const BasicParser & m_owner;
            size_t m_updateCountAtLastRecalc;

            int m_optionIndex = -1;
            StringType m_optionName;
            std::optional<StringType> m_optionArgument;

            int m_positionalIndex = -1;
            std::vector<unsigned> m_positionalSizes;
            
            ValidationData m_validationData;
        };

    private:
        std::vector<Option> m_options;
        std::vector<Positional> m_positionals;
        Tokenizer m_tokenizer;
        std::vector<std::pair<ValidatorFunction, StringType>> m_validators;
        size_t m_updateCount = 0;
        SubCommandMark m_subCommandMark;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(Parser)

}

#endif 
#ifndef HEADER_ARGUM_TYPE_PARSERS_H_INCLUDED
#define HEADER_ARGUM_TYPE_PARSERS_H_INCLUDED



namespace Argum {

    ARGUM_MOD_EXPORTED
    template<class T, StringLike S> 
    requires(std::is_integral_v<T>)
    auto parseIntegral(S && str, int base = 0) -> ARGUM_EXPECTED(CharTypeOf<S>, T) {

        using Char = CharTypeOf<S>;
        using ValidationError = typename BasicParser<Char>::ValidationError;
        using CharConstants = CharConstants<Char>;
        using Messages = Messages<Char>;

        std::basic_string<Char> value(std::forward<S>(str));

        if (value.empty())
            ARGUM_THROW(ValidationError, format(Messages::notANumber(), value));

        T ret;
        Char * endPtr;

        errno = 0;
        if constexpr (std::is_signed_v<T>) {

            if constexpr (sizeof(T) <= sizeof(long)) {
                auto res = CharConstants::toLong(value.data(), &endPtr, base);
                if constexpr (sizeof(T) < sizeof(res)) {
                    if (res < (decltype(res))std::numeric_limits<T>::min() || res > (decltype(res))std::numeric_limits<T>::max())
                        ARGUM_THROW(ValidationError, format(Messages::outOfRange(), value));
                }
                ret = T(res);
            } else  {
                auto res = CharConstants::toLongLong(value.data(), &endPtr, base);
                if constexpr (sizeof(T) < sizeof(res)) {
                    if (res < (decltype(res))std::numeric_limits<T>::min() || res > (decltype(res))std::numeric_limits<T>::max())
                        ARGUM_THROW(ValidationError, format(Messages::outOfRange(), value));
                }
                ret = T(res);
            }

        } else {

            if constexpr (sizeof(T) <= sizeof(unsigned long)) {
                auto res = CharConstants::toULong(value.data(), &endPtr, base);
                if constexpr (sizeof(T) < sizeof(res)) {
                    if (res > (decltype(res))std::numeric_limits<T>::max())
                        ARGUM_THROW(ValidationError, format(Messages::outOfRange(), value));
                }
                ret = T(res);
            } else  {
                auto res = CharConstants::toULongLong(value.data(), &endPtr, base);
                if constexpr (sizeof(T) < sizeof(res)) {
                    if (res > (decltype(res))std::numeric_limits<T>::max())
                        ARGUM_THROW(ValidationError, format(Messages::outOfRange(), value));
                }
                ret = T(res);
            }
        }
        if (errno == ERANGE)
            ARGUM_THROW(ValidationError, format(Messages::outOfRange(), value));
        for ( ; endPtr != value.data() + value.size(); ++endPtr) {
            if (!CharConstants::isSpace(*endPtr))
                ARGUM_THROW(ValidationError, format(Messages::notANumber(), value));
        }
        
        return ret;
    }

    ARGUM_MOD_EXPORTED
    template<class T, StringLike S> 
    requires(std::is_floating_point_v<T>)
    auto parseFloatingPoint(S && str) -> ARGUM_EXPECTED(CharTypeOf<S>, T) {

        using Char = CharTypeOf<S>;
        using ValidationError = typename BasicParser<Char>::ValidationError;
        using CharConstants = CharConstants<Char>;
        using Messages = Messages<Char>;

        std::basic_string<Char> value(std::forward<S>(str));

        if (value.empty())
            ARGUM_THROW(ValidationError, format(Messages::notANumber(), value));

        T ret;
        Char * endPtr;

        errno = 0;
        if constexpr (std::is_same_v<T, float>) {
            ret = CharConstants::toFloat(value.data(), &endPtr);
        } else if constexpr (std::is_same_v<T, double>) {
            ret = CharConstants::toDouble(value.data(), &endPtr);
        } else if constexpr (std::is_same_v<T, long double>) {
            ret = CharConstants::toLongDouble(value.data(), &endPtr);
        }

        if (errno == ERANGE)
            ARGUM_THROW(ValidationError, format(Messages::outOfRange(), value));
        for ( ; endPtr != value.data() + value.size(); ++endPtr) {
            if (!CharConstants::isSpace(*endPtr))
                ARGUM_THROW(ValidationError, format(Messages::notANumber(), value));
        }
        
        return ret;
    }

    ARGUM_MOD_EXPORTED
    template<Character Char>
    class BasicChoiceParser {
    public:
        using StringViewType = std::basic_string_view<Char>;

        struct Settings {
            bool caseSensitive = false;
            bool allowElse = false;
        };
    private:
        using StringType = std::basic_string<Char>;
        using RegexType = std::basic_regex<Char>;
        using CharConstants = Argum::CharConstants<Char>;
        using Messages = Argum::Messages<Char>;
        using ValidationError = typename BasicParser<Char>::ValidationError;
    public:
        BasicChoiceParser(Settings settings = Settings()) : m_options(std::regex_constants::ECMAScript) {
            if (!settings.caseSensitive)
                m_options |= std::regex_constants::icase;
            this->m_allowElse = settings.allowElse;
        }

        template<StringLikeOf<Char> First, StringLikeOf<Char>... Rest>
        auto addChoice(First && first, Rest && ...rest) {
            StringType exp = this->escape(std::forward<First>(first));
            if (!this->m_description.empty())
                this->m_description += Messages::listJoiner();
            this->m_description += std::forward<First>(first);
            this->combine(exp, this->m_description, std::forward<Rest>(rest)...);
            this->m_choices.emplace_back(exp, m_options);
        }

        auto addChoice(std::initializer_list<const Char *> values) {
            auto it = values.begin();
            if (it == values.end())
                ARGUM_INVALID_ARGUMENT("choices list cannot be empty");
            StringType exp = this->escape(*it);
            if (!this->m_description.empty())
                this->m_description += Messages::listJoiner();
            this->m_description += *it;
            for(++it; it != values.end(); ++it)
                this->combine(exp, this->m_description, *it);
            this->m_choices.emplace_back(exp, m_options);
        }

        auto parse(StringViewType value) const -> ARGUM_EXPECTED(Char, size_t) {

            auto it = std::find_if(this->m_choices.begin(), this->m_choices.end(), [&](const RegexType & regex) {
                return regex_match(value.begin(), value.end(), regex);
            });
            if (it == this->m_choices.end()) {
                if (!m_allowElse)
                    ARGUM_THROW(ValidationError, format(Messages::notAValidChoice(), value, this->m_description));
                return this->m_choices.size();
            }
            return it - this->m_choices.begin();
        }

        auto description() const -> StringType {
            return CharConstants::braceOpen + this->m_description + CharConstants::braceClose;
        }
    private:
        auto combine(StringType &, StringType &) {
            return;
        }

        template<StringLikeOf<Char> First, StringLikeOf<Char>... Rest>
        auto combine(StringType & exp, StringType & desc, First && first, Rest && ...rest) {
            exp += CharConstants::pipe;
            exp += this->escape(std::forward<First>(first));
            desc += Messages::listJoiner();
            desc += std::forward<First>(first);
            this->combine(exp, desc, std::forward<Rest>(rest)...);
        }

        template<StringLikeOf<Char> Arg>
        static auto escape(Arg && arg) -> StringType {
            
            StringViewType view(arg);
            if (view.empty())
                ARGUM_INVALID_ARGUMENT("choice cannot be empty");
            
            StringType ret;
            regex_replace(std::back_inserter(ret), view.begin(), view.end(), 
                          BasicChoiceParser::s_needToBeEscaped, BasicChoiceParser::s_escaped,
                          std::regex_constants::match_default | std::regex_constants::format_sed);
            return ret;
        }
    private:
        std::vector<RegexType> m_choices;
        StringType m_description;
        std::regex_constants::syntax_option_type m_options;
        bool m_allowElse;

        static inline const RegexType s_needToBeEscaped{CharConstants::mustEscapeInRegex};
        static inline const StringType s_escaped = CharConstants::regexEscapeReplacement;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(ChoiceParser)

    ARGUM_MOD_EXPORTED
    template<Character Char>
    class BasicBooleanParser {
    public:
        BasicBooleanParser() {
            this->m_impl.addChoice(CharConstants<Char>::falseNames);
            this->m_impl.addChoice(CharConstants<Char>::trueNames);
        }

        auto parse(std::basic_string_view<Char> value) const -> ARGUM_EXPECTED(Char, bool) {
            ARGUM_CHECK_RESULT(auto ret, this->m_impl.parse(value));
            return bool(ret);
        }
    private:
        BasicChoiceParser<Char> m_impl;
    };

    ARGUM_DECLARE_FRIENDLY_NAMES(BooleanParser)

}

#endif
#ifndef HEADER_ARGUM_DETECT_COLOR_H_INCLUDED
#define HEADER_ARGUM_DETECT_COLOR_H_INCLUDED



#if !defined(_WIN32) && __has_include(<unistd.h>)
    #define ARGUM_HAS_UNISTD_H
#endif

#if !defined(_WIN32) 
    #if __has_include(<sys/ioctl.h>)
    #endif
    #if __has_include(<termios.h>)
    #endif
    #ifdef TIOCGWINSZ
        #define ARGUM_HAS_TIOCGWINSZ
    #endif
#endif

#ifdef _WIN32
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #ifndef WIN32_LEAN_AND_MEAN
        #define WIN32_LEAN_AND_MEAN
    #endif
    #ifdef min
        #undef min
    #endif
    #ifdef max
        #undef max
    #endif
#endif

namespace Argum {
    
    ARGUM_MOD_EXPORTED
    enum class ColorStatus {
        unknown,
        forbidden,
        allowed,    
        required    
    };

    // Detect ColorStatus from the environment
    // Logic taken from combination of:
    // https://no-color.org
    // https://force-color.org
    // https://bixense.com/clicolors/
    // https://github.com/chalk/supports-color/blob/main/index.js
    // https://gist.github.com/scop/4d5902b98f0503abec3fcbb00b38aec3
    // https://andrey-zherikov.github.io/argparse/ansi-coloring-and-styling.html#heuristic
    ARGUM_MOD_EXPORTED
    inline auto environmentColorStatus() -> ColorStatus {
        using namespace std;
        using namespace std::literals;

        if (auto val = getenv("NO_COLOR"); val && *val)
            return ColorStatus::forbidden;

        if (auto val = getenv("FORCE_COLOR"); val && *val)
            return ColorStatus::required;

        if (auto val = getenv("CLICOLOR_FORCE"); val && *val) {
            if (strcmp(val, "0") == 0 || strcmp(val, "false") == 0)
                return ColorStatus::forbidden;
            return ColorStatus::required;
        }

        if (auto val = getenv("CLICOLOR"); val && *val) {
            if (strcmp(val, "0") == 0 || strcmp(val, "false") == 0)
                return ColorStatus::forbidden;
            return ColorStatus::allowed;
        }

#ifdef _WIN32
        if (auto val = getenv("WT_SESSION"); val && *val)
            return ColorStatus::allowed;
#endif

        if (auto val = getenv("COLORTERM"); val && *val)
            return ColorStatus::allowed;

        if (auto val = getenv("TERM")) {
            string_view term(val);

            if (term == "dumb")
                return ColorStatus::forbidden;

            for (auto & exact: {"wezterm"sv}) {
                if (term == exact)
                    return ColorStatus::allowed;
            }
            for (auto & start: {"screen"sv, "xterm"sv, "vt100"sv, "vt220"sv, "rxvt"sv,
                                "gnome"sv, "konsole"sv, "kterm"sv, "alacritty"sv, "console"sv}) {
                if (term.size() >= start.size() && term.substr(0, start.size()) == start)
                    return ColorStatus::allowed;
            }
            for (auto & inside: {"color"sv, "ansi"sv, "cygwin"sv, "linux"sv}) {
                if (term.find(inside) != term.npos)
                    return ColorStatus::allowed;
            }
            for (auto & end: {"-256"sv}) {
                if (term.size() >= end.size() && term.substr(term.size() - end.size()) == end)
                    return ColorStatus::allowed;
            }
        }

        return ColorStatus::unknown;
    }

    ARGUM_MOD_EXPORTED
    inline bool shouldUseColor(ColorStatus envColorStatus, FILE * fp) {
        if (envColorStatus == ColorStatus::required)
            return true;
        if (envColorStatus == ColorStatus::forbidden)
            return false;

#if defined(ARGUM_HAS_UNISTD_H)
        if (envColorStatus == ColorStatus::unknown)
            return false;

        return isatty(fileno(fp));

#elif defined(_WIN32)
        int desc = _fileno(fp);
        if (desc < 0)
            return false;

        if (!_isatty(desc))
            return false;

        if (envColorStatus == ColorStatus::unknown) {
        
            HANDLE h = HANDLE(_get_osfhandle(_fileno(fp)));
            if (h == INVALID_HANDLE_VALUE)
                return false;

            DWORD dwMode = 0;
            if (!GetConsoleMode(h, &dwMode))
                return false;

            if (!(dwMode & ENABLE_VIRTUAL_TERMINAL_PROCESSING))
                return false;
        }

        return true;
#else
        return false;
#endif
    }

    ARGUM_MOD_EXPORTED
    inline auto colorizerForFile(ColorStatus envColorStatus, FILE * fp) -> Colorizer {
        if (shouldUseColor(envColorStatus, fp))
            return defaultColorizer();
        return {};
    }

    ARGUM_MOD_EXPORTED
    inline auto wideColorizerForFile(ColorStatus envColorStatus, FILE * fp) -> WColorizer {
        if (shouldUseColor(envColorStatus, fp))
            return defaultWColorizer();
        return {};
    }

    ARGUM_MOD_EXPORTED
    inline unsigned terminalWidth(FILE * fp) {
        unsigned fallback = std::numeric_limits<unsigned>::max();

#if defined(ARGUM_HAS_UNISTD_H) && defined(ARGUM_HAS_TIOCGWINSZ) 
        struct winsize ws{};

        int desc = fileno(fp);
        if (isatty(desc) &&
            ioctl(desc, TIOCGWINSZ, &ws) == 0 &&
            ws.ws_col > 0) {

            return ws.ws_col;
        }

        if (auto * cols = getenv("COLUMNS")) {
            char * end;
            auto val = CharConstants<char>::toULong(cols, &end, 10);
            if (val > 0 && val < std::numeric_limits<unsigned>::max() &&
                end == cols + strlen(cols)) {

                return unsigned(val);
            }
        }

        return fallback;

#elif defined(_WIN32)

        int desc = _fileno(fp);
        if (desc < 0)
            return fallback;

        if (!_isatty(desc))
            return fallback;

        HANDLE h = HANDLE(_get_osfhandle(desc));
        if (h == INVALID_HANDLE_VALUE)
            return fallback;

        CONSOLE_SCREEN_BUFFER_INFO csbi;
        if (!GetConsoleScreenBufferInfo(h, &csbi))
            return fallback;

        return csbi.srWindow.Right - csbi.srWindow.Left + 1;

#else
        return fallback;
#endif
    }

}

#endif

#endif


#endif 

