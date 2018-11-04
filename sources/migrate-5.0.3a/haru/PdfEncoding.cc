/*
 * << H a r u --free pdf library >> -- PdfEncodingDef.h
 *
 * Copyright (c) 1999-2003 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 */

#include "libharu.h"
#include <new>

static const unsigned char UNICODE_HEADER[PDF_UNICODE_HEADER_LEN] = {
    0xFE, 0xFF
};

struct pdf_unicode_gryph_pair_struct {
    unsigned short unicode;
    const char* gryph_name;
};
typedef pdf_unicode_gryph_pair_struct pdf_unicode_gryph_pair;

/*----------------------------------------------------------------------------*/

struct pdf_encode_name_pair {
    pdf_encoding encode;
    char* encode_name;
};

static const pdf_encode_name_pair pdf_encoding_names[5] = {
  {PDF_STANDARD_ENCODING, (char *) "StandardEncoding"},
  {PDF_MAC_ROMAN_ENCODING, (char *)"MacRomanEncoding"},
  {PDF_WIN_ANSI_ENCODING, (char *) "WinAnsiEncoding"},
  {PDF_FONT_SPECIFIC, (char *) "FontSpecific"}
};

/*----------------------------------------------------------------------------*/

PdfEncodingDef::PdfEncodingDef()
    : PdfAutoPtrObject()
{
    fFirstChar = 0;
    fLastChar = 0;
    fBaseEncoding = PDF_ENCODING_EOF;
    fDifferences = NULL;

    for (int i = 0; i < 255; i++)
        fUnicodeArray[i] = 0x0000;
        
    PDF_DEBUG_PRINT(("++ [%x] PdfEncodingDef new.\n", (int)this));
}

void
PdfEncodingDef::AddDifference(unsigned int code, const char* char_name)
{
    if (fDifferences == NULL) {
        fDifferences = new PdfList();
        PDF_DEBUG_PRINT(("++ [%x] fDifferences new.\n", 
                    (int)fDifferences));
    }
        
    pdf_differences_element* diff = new pdf_differences_element;
    PDF_DEBUG_PRINT(("++ [%x] pdf_differences_element new.\n", (int)diff));
    diff->code = code;
    diff->char_name = char_name;

    PDF_DEBUG_PRINT(("PdfEncodingDef::AddDifference [%d], [%d]\n",
                code, (int)char_name));
    if (fDifferences->AddItem(diff) != true) {
        delete diff;
        throw PdfException(PDF_ERR_MALLOC, 
                "ERROR: PdfEncodingDef::AddDifference");
    }
}

int
PdfEncodingDef::GetCharCode(unsigned short unicode)
{
    for (unsigned int i = 0; i <= fLastChar - fFirstChar; i++) {
        if (fUnicodeArray[i] == unicode)
            return i + fFirstChar;
    }
    return -1;
}

int
PdfEncodingDef::GetCharCode(const char* char_name)
{
    unsigned int unicode = GryphNameToUnicode(char_name);
    return GetCharCode(unicode);
}

const char*
PdfEncodingDef::GetCharName(unsigned int code)
{
    if (fFirstChar > code || code > fLastChar)
        return PDF_CHAR_NOTDEF;
    
    unsigned int unicode = fUnicodeArray[code - fFirstChar];
    return UnicodeToGryphName(unicode);
}

unsigned short
PdfEncodingDef::GetUnicode(unsigned int code)
{
    return ((fFirstChar <= code && code <= fLastChar) ?
       fUnicodeArray[code - fFirstChar] : 0x0000);
}

PdfEncodingDef::~PdfEncodingDef()
{
    ClearDefferences();
    
    PDF_DEBUG_PRINT(("++ [%x] fDifferences delete.\n", (int)fDifferences));
    delete fDifferences;

    PDF_DEBUG_PRINT(("++ [%x] PdfEncodingDef delete.\n", (int)this));
}

void
PdfEncodingDef::ClearDefferences()
{
    if (fDifferences != NULL) {
        for (int i = 0; i < fDifferences->CountItems(); i++) {
            PDF_DEBUG_PRINT(("++ [%x] Differences delete.\n",
                        (int)fDifferences->ItemAt(i)));
            delete (pdf_differences_element*)fDifferences->ItemAt(i);
        }
    }
}

PdfObject*
PdfEncodingDef::GetEncoding(PdfXref *xref)
{
    if (fDifferences == NULL || fDifferences->CountItems() == 0) {
        PdfName* name = new PdfName(PdfEncodingToString(fBaseEncoding));
        return name;
    } else {

        /* create dictionary for encoding. */
        PdfDictionary *dict = new PdfDictionary(xref);
        try {
            dict->AddElement("Type", new PdfName("Encoding"));

            if (fBaseEncoding != PDF_FONT_SPECIFIC) {
                const char* base_name = PdfEncodingToString(fBaseEncoding);

                dict->AddElement("BaseEncoding", new PdfName(base_name));
            }
            PdfArray* array = new PdfArray(xref);

            try {
                for (int i = 0; i < fDifferences->CountItems(); i++) {
                    pdf_differences_element* element =
                        (pdf_differences_element*)fDifferences->ItemAt(i);
                    array->Add(new PdfNumber(element->code));
                    array->Add(new PdfName(element->char_name));
                }
            } catch (...) {
                delete array;
                throw;
            }
            dict->AddElement("Differences", array);
        } catch (...) {
            delete dict;
            throw;
        }
        return dict;
    }
}

int
PdfEncodingDef::ToUnicode(const char* src, unsigned char* dst, int* len)
{
    int srclen = strlen(src);
    int dstlen;
    int j = 0;

    if (dst != NULL && len != NULL)
        dstlen = *len - 1;
    else
        dstlen = -1;

    if (dstlen > 2) {
        memcpy(dst, UNICODE_HEADER, PDF_UNICODE_HEADER_LEN);
        j = 2;
    } else 
        return -1;

    const unsigned char* usrc = (const unsigned char*)src;
    for (int i = 0; i < srclen; i++) {
        unsigned int unicode;
        unsigned char tmp[2];

        unicode = GetUnicode(usrc[i]);
        tmp[0] = unicode >> 8;
        tmp[1] = unicode;

        if (j < dstlen) {
            memcpy(dst + j, tmp, 2);
            j += 2;
        }
    }
    *len = j;
    return j;
}

void
PdfEncodingDef::SetUnicodeArray(unsigned int first_char,
        unsigned int last_char, const unsigned short* array)
{
    fFirstChar = first_char;
    fLastChar = last_char;

    ClearDefferences();
    
    for (unsigned int i = 0; i <= fLastChar - fFirstChar; i++)
        fUnicodeArray[i] = array[i];
}

void
PdfEncodingDef::OverrideUnicodeArray(const unsigned short* array)
{
    ClearDefferences();

    for (unsigned int i = 0; i <= fLastChar - fFirstChar; i++) {
        if (fUnicodeArray[i] != array[i]) {
            fUnicodeArray[i] = array[i];
            AddDifference(i + fFirstChar, UnicodeToGryphName(array[i]));
        }
    }
}

/*---------------------------------------------------------------------------*/
/*------ utility routine ----------------------------------------------------*/

static const pdf_unicode_gryph_pair UNICODE_GRYPH_NAME_MAP[] = {
    {0x0000, PDF_CHAR_NOTDEF},
    {0x0020, "space"},
    {0x0021, "exclam"},
    {0x0022, "quotedbl"},
    {0x0023, "numbersign"},
    {0x0024, "dollar"},
    {0x0025, "percent"},
    {0x0026, "ampersand"},
    {0x0027, "quotesingle"},
    {0x0028, "parenleft"},
    {0x0029, "parenright"},
    {0x002A, "asterisk"},
    {0x002B, "plus"},
    {0x002C, "comma"},
    {0x002D, "hyphen"},
    {0x002E, "period"},
    {0x002F, "slash"},
    {0x0030, "zero"},
    {0x0031, "one"},
    {0x0032, "two"},
    {0x0033, "three"},
    {0x0034, "four"},
    {0x0035, "five"},
    {0x0036, "six"},
    {0x0037, "seven"},
    {0x0038, "eight"},
    {0x0039, "nine"},
    {0x003A, "colon"},
    {0x003B, "semicolon"},
    {0x003C, "less"},
    {0x003D, "equal"},
    {0x003E, "greater"},
    {0x003F, "question"},
    {0x0040, "at"},
    {0x0041, "A"},
    {0x0042, "B"},
    {0x0043, "C"},
    {0x0044, "D"},
    {0x0045, "E"},
    {0x0046, "F"},
    {0x0047, "G"},
    {0x0048, "H"},
    {0x0049, "I"},
    {0x004A, "J"},
    {0x004B, "K"},
    {0x004C, "L"},
    {0x004D, "M"},
    {0x004E, "N"},
    {0x004F, "O"},
    {0x0050, "P"},
    {0x0051, "Q"},
    {0x0052, "R"},
    {0x0053, "S"},
    {0x0054, "T"},
    {0x0055, "U"},
    {0x0056, "V"},
    {0x0057, "W"},
    {0x0058, "X"},
    {0x0059, "Y"},
    {0x005A, "Z"},
    {0x005B, "bracketleft"},
    {0x005C, "backslash"},
    {0x005D, "bracketright"},
    {0x005E, "asciicircum"},
    {0x005F, "underscore"},
    {0x0060, "grave"},
    {0x0061, "a"},
    {0x0062, "b"},
    {0x0063, "c"},
    {0x0064, "d"},
    {0x0065, "e"},
    {0x0066, "f"},
    {0x0067, "g"},
    {0x0068, "h"},
    {0x0069, "i"},
    {0x006A, "j"},
    {0x006B, "k"},
    {0x006C, "l"},
    {0x006D, "m"},
    {0x006E, "n"},
    {0x006F, "o"},
    {0x0070, "p"},
    {0x0071, "q"},
    {0x0072, "r"},
    {0x0073, "s"},
    {0x0074, "t"},
    {0x0075, "u"},
    {0x0076, "v"},
    {0x0077, "w"},
    {0x0078, "x"},
    {0x0079, "y"},
    {0x007A, "z"},
    {0x007B, "braceleft"},
    {0x007C, "bar"},
    {0x007D, "braceright"},
    {0x007E, "asciitilde"},
    {0x00A0, "space"},
    {0x00A1, "exclamdown"},
    {0x00A2, "cent"},
    {0x00A3, "sterling"},
    {0x00A4, "currency"},
    {0x00A5, "yen"},
    {0x00A6, "brokenbar"},
    {0x00A7, "section"},
    {0x00A8, "dieresis"},
    {0x00A9, "copyright"},
    {0x00AA, "ordfeminine"},
    {0x00AB, "guillemotleft"},
    {0x00AC, "logicalnot"},
    {0x00AD, "hyphen"},
    {0x00AE, "registered"},
    {0x00AF, "macron"},
    {0x00B0, "degree"},
    {0x00B1, "plusminus"},
    {0x00B2, "twosuperior"},
    {0x00B3, "threesuperior"},
    {0x00B4, "acute"},
    {0x00B5, "mu"},
    {0x00B6, "paragraph"},
    {0x00B7, "periodcentered"},
    {0x00B8, "cedilla"},
    {0x00B9, "onesuperior"},
    {0x00BA, "ordmasculine"},
    {0x00BB, "guillemotright"},
    {0x00BC, "onequarter"},
    {0x00BD, "onehalf"},
    {0x00BE, "threequarters"},
    {0x00BF, "questiondown"},
    {0x00C0, "Agrave"},
    {0x00C1, "Aacute"},
    {0x00C2, "Acircumflex"},
    {0x00C3, "Atilde"},
    {0x00C4, "Adieresis"},
    {0x00C5, "Aring"},
    {0x00C6, "AE"},
    {0x00C7, "Ccedilla"},
    {0x00C8, "Egrave"},
    {0x00C9, "Eacute"},
    {0x00CA, "Ecircumflex"},
    {0x00CB, "Edieresis"},
    {0x00CC, "Igrave"},
    {0x00CD, "Iacute"},
    {0x00CE, "Icircumflex"},
    {0x00CF, "Idieresis"},
    {0x00D0, "Eth"},
    {0x00D1, "Ntilde"},
    {0x00D2, "Ograve"},
    {0x00D3, "Oacute"},
    {0x00D4, "Ocircumflex"},
    {0x00D5, "Otilde"},
    {0x00D6, "Odieresis"},
    {0x00D7, "multiply"},
    {0x00D8, "Oslash"},
    {0x00D9, "Ugrave"},
    {0x00DA, "Uacute"},
    {0x00DB, "Ucircumflex"},
    {0x00DC, "Udieresis"},
    {0x00DD, "Yacute"},
    {0x00DE, "Thorn"},
    {0x00DF, "germandbls"},
    {0x00E0, "agrave"},
    {0x00E1, "aacute"},
    {0x00E2, "acircumflex"},
    {0x00E3, "atilde"},
    {0x00E4, "adieresis"},
    {0x00E5, "aring"},
    {0x00E6, "ae"},
    {0x00E7, "ccedilla"},
    {0x00E8, "egrave"},
    {0x00E9, "eacute"},
    {0x00EA, "ecircumflex"},
    {0x00EB, "edieresis"},
    {0x00EC, "igrave"},
    {0x00ED, "iacute"},
    {0x00EE, "icircumflex"},
    {0x00EF, "idieresis"},
    {0x00F0, "eth"},
    {0x00F1, "ntilde"},
    {0x00F2, "ograve"},
    {0x00F3, "oacute"},
    {0x00F4, "ocircumflex"},
    {0x00F5, "otilde"},
    {0x00F6, "odieresis"},
    {0x00F7, "divide"},
    {0x00F8, "oslash"},
    {0x00F9, "ugrave"},
    {0x00FA, "uacute"},
    {0x00FB, "ucircumflex"},
    {0x00FC, "udieresis"},
    {0x00FD, "yacute"},
    {0x00FE, "thorn"},
    {0x00FF, "ydieresis"},
    {0x0100, "Amacron"},
    {0x0101, "amacron"},
    {0x0102, "Abreve"},
    {0x0103, "abreve"},
    {0x0104, "Aogonek"},
    {0x0105, "aogonek"},
    {0x0106, "Cacute"},
    {0x0107, "cacute"},
    {0x0108, "Ccircumflex"},
    {0x0109, "ccircumflex"},
    {0x010A, "Cdotaccent"},
    {0x010B, "cdotaccent"},
    {0x010C, "Ccaron"},
    {0x010D, "ccaron"},
    {0x010E, "Dcaron"},
    {0x010F, "dcaron"},
    {0x0110, "Dcroat"},
    {0x0111, "dcroat"},
    {0x0112, "Emacron"},
    {0x0113, "emacron"},
    {0x0114, "Ebreve"},
    {0x0115, "ebreve"},
    {0x0116, "Edotaccent"},
    {0x0117, "edotaccent"},
    {0x0118, "Eogonek"},
    {0x0119, "eogonek"},
    {0x011A, "Ecaron"},
    {0x011B, "ecaron"},
    {0x011C, "Gcircumflex"},
    {0x011D, "gcircumflex"},
    {0x011E, "Gbreve"},
    {0x011F, "gbreve"},
    {0x0120, "Gdotaccent"},
    {0x0121, "gdotaccent"},
    {0x0122, "Gcommaaccent"},
    {0x0123, "gcommaaccent"},
    {0x0124, "Hcircumflex"},
    {0x0125, "hcircumflex"},
    {0x0126, "Hbar"},
    {0x0127, "hbar"},
    {0x0128, "Itilde"},
    {0x0129, "itilde"},
    {0x012A, "Imacron"},
    {0x012B, "imacron"},
    {0x012C, "Ibreve"},
    {0x012D, "ibreve"},
    {0x012E, "Iogonek"},
    {0x012F, "iogonek"},
    {0x0130, "Idotaccent"},
    {0x0131, "dotlessi"},
    {0x0132, "IJ"},
    {0x0133, "ij"},
    {0x0134, "Jcircumflex"},
    {0x0135, "jcircumflex"},
    {0x0136, "Kcommaaccent"},
    {0x0137, "kcommaaccent"},
    {0x0138, "kgreenlandic"},
    {0x0139, "Lacute"},
    {0x013A, "lacute"},
    {0x013B, "Lcommaaccent"},
    {0x013C, "lcommaaccent"},
    {0x013D, "Lcaron"},
    {0x013E, "lcaron"},
    {0x013F, "Ldot"},
    {0x0140, "ldot"},
    {0x0141, "Lslash"},
    {0x0142, "lslash"},
    {0x0143, "Nacute"},
    {0x0144, "nacute"},
    {0x0145, "Ncommaaccent"},
    {0x0146, "ncommaaccent"},
    {0x0147, "Ncaron"},
    {0x0148, "ncaron"},
    {0x0149, "napostrophe"},
    {0x014A, "Eng"},
    {0x014B, "eng"},
    {0x014C, "Omacron"},
    {0x014D, "omacron"},
    {0x014E, "Obreve"},
    {0x014F, "obreve"},
    {0x0150, "Ohungarumlaut"},
    {0x0151, "ohungarumlaut"},
    {0x0152, "OE"},
    {0x0153, "oe"},
    {0x0154, "Racute"},
    {0x0155, "racute"},
    {0x0156, "Rcommaaccent"},
    {0x0157, "rcommaaccent"},
    {0x0158, "Rcaron"},
    {0x0159, "rcaron"},
    {0x015A, "Sacute"},
    {0x015B, "sacute"},
    {0x015C, "Scircumflex"},
    {0x015D, "scircumflex"},
    {0x015E, "Scedilla"},
    {0x015F, "scedilla"},
    {0x0160, "Scaron"},
    {0x0161, "scaron"},
    {0x0162, "Tcommaaccent"},
    {0x0163, "tcommaaccent"},
    {0x0164, "Tcaron"},
    {0x0165, "tcaron"},
    {0x0166, "Tbar"},
    {0x0167, "tbar"},
    {0x0168, "Utilde"},
    {0x0169, "utilde"},
    {0x016A, "Umacron"},
    {0x016B, "umacron"},
    {0x016C, "Ubreve"},
    {0x016D, "ubreve"},
    {0x016E, "Uring"},
    {0x016F, "uring"},
    {0x0170, "Uhungarumlaut"},
    {0x0171, "uhungarumlaut"},
    {0x0172, "Uogonek"},
    {0x0173, "uogonek"},
    {0x0174, "Wcircumflex"},
    {0x0175, "wcircumflex"},
    {0x0176, "Ycircumflex"},
    {0x0177, "ycircumflex"},
    {0x0178, "Ydieresis"},
    {0x0179, "Zacute"},
    {0x017A, "zacute"},
    {0x017B, "Zdotaccent"},
    {0x017C, "zdotaccent"},
    {0x017D, "Zcaron"},
    {0x017E, "zcaron"},
    {0x017F, "longs"},
    {0x0192, "florin"},
    {0x01A0, "Ohorn"},
    {0x01A1, "ohorn"},
    {0x01AF, "Uhorn"},
    {0x01B0, "uhorn"},
    {0x01E6, "Gcaron"},
    {0x01E7, "gcaron"},
    {0x01FA, "Aringacute"},
    {0x01FB, "aringacute"},
    {0x01FC, "AEacute"},
    {0x01FD, "aeacute"},
    {0x01FE, "Oslashacute"},
    {0x01FF, "oslashacute"},
    {0x0218, "Scommaaccent"},
    {0x0219, "scommaaccent"},
    {0x021A, "Tcommaaccent"},
    {0x021B, "tcommaaccent"},
    {0x02BC, "afii57929"},
    {0x02BD, "afii64937"},
    {0x02C6, "circumflex"},
    {0x02C7, "caron"},
    {0x02C9, "macron"},
    {0x02D8, "breve"},
    {0x02D9, "dotaccent"},
    {0x02DA, "ring"},
    {0x02DB, "ogonek"},
    {0x02DC, "tilde"},
    {0x02DD, "hungarumlaut"},
    {0x0300, "gravecomb"},
    {0x0301, "acutecomb"},
    {0x0303, "tildecomb"},
    {0x0309, "hookabovecomb"},
    {0x0323, "dotbelowcomb"},
    {0x0384, "tonos"},
    {0x0385, "dieresistonos"},
    {0x0386, "Alphatonos"},
    {0x0387, "anoteleia"},
    {0x0388, "Epsilontonos"},
    {0x0389, "Etatonos"},
    {0x038A, "Iotatonos"},
    {0x038C, "Omicrontonos"},
    {0x038E, "Upsilontonos"},
    {0x038F, "Omegatonos"},
    {0x0390, "iotadieresistonos"},
    {0x0391, "Alpha"},
    {0x0392, "Beta"},
    {0x0393, "Gamma"},
    {0x0394, "Delta"},
    {0x0395, "Epsilon"},
    {0x0396, "Zeta"},
    {0x0397, "Eta"},
    {0x0398, "Theta"},
    {0x0399, "Iota"},
    {0x039A, "Kappa"},
    {0x039B, "Lambda"},
    {0x039C, "Mu"},
    {0x039D, "Nu"},
    {0x039E, "Xi"},
    {0x039F, "Omicron"},
    {0x03A0, "Pi"},
    {0x03A1, "Rho"},
    {0x03A3, "Sigma"},
    {0x03A4, "Tau"},
    {0x03A5, "Upsilon"},
    {0x03A6, "Phi"},
    {0x03A7, "Chi"},
    {0x03A8, "Psi"},
    {0x03A9, "Omega"},
    {0x03AA, "Iotadieresis"},
    {0x03AB, "Upsilondieresis"},
    {0x03AC, "alphatonos"},
    {0x03AD, "epsilontonos"},
    {0x03AE, "etatonos"},
    {0x03AF, "iotatonos"},
    {0x03B0, "upsilondieresistonos"},
    {0x03B1, "alpha"},
    {0x03B2, "beta"},
    {0x03B3, "gamma"},
    {0x03B4, "delta"},
    {0x03B5, "epsilon"},
    {0x03B6, "zeta"},
    {0x03B7, "eta"},
    {0x03B8, "theta"},
    {0x03B9, "iota"},
    {0x03BA, "kappa"},
    {0x03BB, "lambda"},
    {0x03BC, "mu"},
    {0x03BD, "nu"},
    {0x03BE, "xi"},
    {0x03BF, "omicron"},
    {0x03C0, "pi"},
    {0x03C1, "rho"},
    {0x03C2, "sigma1"},
    {0x03C3, "sigma"},
    {0x03C4, "tau"},
    {0x03C5, "upsilon"},
    {0x03C6, "phi"},
    {0x03C7, "chi"},
    {0x03C8, "psi"},
    {0x03C9, "omega"},
    {0x03CA, "iotadieresis"},
    {0x03CB, "upsilondieresis"},
    {0x03CC, "omicrontonos"},
    {0x03CD, "upsilontonos"},
    {0x03CE, "omegatonos"},
    {0x03D1, "theta1"},
    {0x03D2, "Upsilon1"},
    {0x03D5, "phi1"},
    {0x03D6, "omega1"},
    {0x0401, "afii10023"},
    {0x0402, "afii10051"},
    {0x0403, "afii10052"},
    {0x0404, "afii10053"},
    {0x0405, "afii10054"},
    {0x0406, "afii10055"},
    {0x0407, "afii10056"},
    {0x0408, "afii10057"},
    {0x0409, "afii10058"},
    {0x040A, "afii10059"},
    {0x040B, "afii10060"},
    {0x040C, "afii10061"},
    {0x040E, "afii10062"},
    {0x040F, "afii10145"},
    {0x0410, "afii10017"},
    {0x0411, "afii10018"},
    {0x0412, "afii10019"},
    {0x0413, "afii10020"},
    {0x0414, "afii10021"},
    {0x0415, "afii10022"},
    {0x0416, "afii10024"},
    {0x0417, "afii10025"},
    {0x0418, "afii10026"},
    {0x0419, "afii10027"},
    {0x041A, "afii10028"},
    {0x041B, "afii10029"},
    {0x041C, "afii10030"},
    {0x041D, "afii10031"},
    {0x041E, "afii10032"},
    {0x041F, "afii10033"},
    {0x0420, "afii10034"},
    {0x0421, "afii10035"},
    {0x0422, "afii10036"},
    {0x0423, "afii10037"},
    {0x0424, "afii10038"},
    {0x0425, "afii10039"},
    {0x0426, "afii10040"},
    {0x0427, "afii10041"},
    {0x0428, "afii10042"},
    {0x0429, "afii10043"},
    {0x042A, "afii10044"},
    {0x042B, "afii10045"},
    {0x042C, "afii10046"},
    {0x042D, "afii10047"},
    {0x042E, "afii10048"},
    {0x042F, "afii10049"},
    {0x0430, "afii10065"},
    {0x0431, "afii10066"},
    {0x0432, "afii10067"},
    {0x0433, "afii10068"},
    {0x0434, "afii10069"},
    {0x0435, "afii10070"},
    {0x0436, "afii10072"},
    {0x0437, "afii10073"},
    {0x0438, "afii10074"},
    {0x0439, "afii10075"},
    {0x043A, "afii10076"},
    {0x043B, "afii10077"},
    {0x043C, "afii10078"},
    {0x043D, "afii10079"},
    {0x043E, "afii10080"},
    {0x043F, "afii10081"},
    {0x0440, "afii10082"},
    {0x0441, "afii10083"},
    {0x0442, "afii10084"},
    {0x0443, "afii10085"},
    {0x0444, "afii10086"},
    {0x0445, "afii10087"},
    {0x0446, "afii10088"},
    {0x0447, "afii10089"},
    {0x0448, "afii10090"},
    {0x0449, "afii10091"},
    {0x044A, "afii10092"},
    {0x044B, "afii10093"},
    {0x044C, "afii10094"},
    {0x044D, "afii10095"},
    {0x044E, "afii10096"},
    {0x044F, "afii10097"},
    {0x0451, "afii10071"},
    {0x0452, "afii10099"},
    {0x0453, "afii10100"},
    {0x0454, "afii10101"},
    {0x0455, "afii10102"},
    {0x0456, "afii10103"},
    {0x0457, "afii10104"},
    {0x0458, "afii10105"},
    {0x0459, "afii10106"},
    {0x045A, "afii10107"},
    {0x045B, "afii10108"},
    {0x045C, "afii10109"},
    {0x045E, "afii10110"},
    {0x045F, "afii10193"},
    {0x0462, "afii10146"},
    {0x0463, "afii10194"},
    {0x0472, "afii10147"},
    {0x0473, "afii10195"},
    {0x0474, "afii10148"},
    {0x0475, "afii10196"},
    {0x0490, "afii10050"},
    {0x0491, "afii10098"},
    {0x04D9, "afii10846"},
    {0x05B0, "afii57799"},
    {0x05B1, "afii57801"},
    {0x05B2, "afii57800"},
    {0x05B3, "afii57802"},
    {0x05B4, "afii57793"},
    {0x05B5, "afii57794"},
    {0x05B6, "afii57795"},
    {0x05B7, "afii57798"},
    {0x05B8, "afii57797"},
    {0x05B9, "afii57806"},
    {0x05BB, "afii57796"},
    {0x05BC, "afii57807"},
    {0x05BD, "afii57839"},
    {0x05BE, "afii57645"},
    {0x05BF, "afii57841"},
    {0x05C0, "afii57842"},
    {0x05C1, "afii57804"},
    {0x05C2, "afii57803"},
    {0x05C3, "afii57658"},
    {0x05D0, "afii57664"},
    {0x05D1, "afii57665"},
    {0x05D2, "afii57666"},
    {0x05D3, "afii57667"},
    {0x05D4, "afii57668"},
    {0x05D5, "afii57669"},
    {0x05D6, "afii57670"},
    {0x05D7, "afii57671"},
    {0x05D8, "afii57672"},
    {0x05D9, "afii57673"},
    {0x05DA, "afii57674"},
    {0x05DB, "afii57675"},
    {0x05DC, "afii57676"},
    {0x05DD, "afii57677"},
    {0x05DE, "afii57678"},
    {0x05DF, "afii57679"},
    {0x05E0, "afii57680"},
    {0x05E1, "afii57681"},
    {0x05E2, "afii57682"},
    {0x05E3, "afii57683"},
    {0x05E4, "afii57684"},
    {0x05E5, "afii57685"},
    {0x05E6, "afii57686"},
    {0x05E7, "afii57687"},
    {0x05E8, "afii57688"},
    {0x05E9, "afii57689"},
    {0x05EA, "afii57690"},
    {0x05F0, "afii57716"},
    {0x05F1, "afii57717"},
    {0x05F2, "afii57718"},
    {0x060C, "afii57388"},
    {0x061B, "afii57403"},
    {0x061F, "afii57407"},
    {0x0621, "afii57409"},
    {0x0622, "afii57410"},
    {0x0623, "afii57411"},
    {0x0624, "afii57412"},
    {0x0625, "afii57413"},
    {0x0626, "afii57414"},
    {0x0627, "afii57415"},
    {0x0628, "afii57416"},
    {0x0629, "afii57417"},
    {0x062A, "afii57418"},
    {0x062B, "afii57419"},
    {0x062C, "afii57420"},
    {0x062D, "afii57421"},
    {0x062E, "afii57422"},
    {0x062F, "afii57423"},
    {0x0630, "afii57424"},
    {0x0631, "afii57425"},
    {0x0632, "afii57426"},
    {0x0633, "afii57427"},
    {0x0634, "afii57428"},
    {0x0635, "afii57429"},
    {0x0636, "afii57430"},
    {0x0637, "afii57431"},
    {0x0638, "afii57432"},
    {0x0639, "afii57433"},
    {0x063A, "afii57434"},
    {0x0640, "afii57440"},
    {0x0641, "afii57441"},
    {0x0642, "afii57442"},
    {0x0643, "afii57443"},
    {0x0644, "afii57444"},
    {0x0645, "afii57445"},
    {0x0646, "afii57446"},
    {0x0647, "afii57470"},
    {0x0648, "afii57448"},
    {0x0649, "afii57449"},
    {0x064A, "afii57450"},
    {0x064B, "afii57451"},
    {0x064C, "afii57452"},
    {0x064D, "afii57453"},
    {0x064E, "afii57454"},
    {0x064F, "afii57455"},
    {0x0650, "afii57456"},
    {0x0651, "afii57457"},
    {0x0652, "afii57458"},
    {0x0660, "afii57392"},
    {0x0661, "afii57393"},
    {0x0662, "afii57394"},
    {0x0663, "afii57395"},
    {0x0664, "afii57396"},
    {0x0665, "afii57397"},
    {0x0666, "afii57398"},
    {0x0667, "afii57399"},
    {0x0668, "afii57400"},
    {0x0669, "afii57401"},
    {0x066A, "afii57381"},
    {0x066D, "afii63167"},
    {0x0679, "afii57511"},
    {0x067E, "afii57506"},
    {0x0686, "afii57507"},
    {0x0688, "afii57512"},
    {0x0691, "afii57513"},
    {0x0698, "afii57508"},
    {0x06A4, "afii57505"},
    {0x06AF, "afii57509"},
    {0x06BA, "afii57514"},
    {0x06D2, "afii57519"},
    {0x06D5, "afii57534"},
    {0x1E80, "Wgrave"},
    {0x1E81, "wgrave"},
    {0x1E82, "Wacute"},
    {0x1E83, "wacute"},
    {0x1E84, "Wdieresis"},
    {0x1E85, "wdieresis"},
    {0x1EF2, "Ygrave"},
    {0x1EF3, "ygrave"},
    {0x200C, "afii61664"},
    {0x200D, "afii301"},
    {0x200E, "afii299"},
    {0x200F, "afii300"},
    {0x2012, "figuredash"},
    {0x2013, "endash"},
    {0x2014, "emdash"},
    {0x2015, "afii00208"},
    {0x2017, "underscoredbl"},
    {0x2018, "quoteleft"},
    {0x2019, "quoteright"},
    {0x201A, "quotesinglbase"},
    {0x201B, "quotereversed"},
    {0x201C, "quotedblleft"},
    {0x201D, "quotedblright"},
    {0x201E, "quotedblbase"},
    {0x2020, "dagger"},
    {0x2021, "daggerdbl"},
    {0x2022, "bullet"},
    {0x2024, "onedotenleader"},
    {0x2025, "twodotenleader"},
    {0x2026, "ellipsis"},
    {0x202C, "afii61573"},
    {0x202D, "afii61574"},
    {0x202E, "afii61575"},
    {0x2030, "perthousand"},
    {0x2032, "minute"},
    {0x2033, "second"},
    {0x2039, "guilsinglleft"},
    {0x203A, "guilsinglright"},
    {0x203C, "exclamdbl"},
    {0x2044, "fraction"},
    {0x2070, "zerosuperior"},
    {0x2074, "foursuperior"},
    {0x2075, "fivesuperior"},
    {0x2076, "sixsuperior"},
    {0x2077, "sevensuperior"},
    {0x2078, "eightsuperior"},
    {0x2079, "ninesuperior"},
    {0x207D, "parenleftsuperior"},
    {0x207E, "parenrightsuperior"},
    {0x207F, "nsuperior"},
    {0x2080, "zeroinferior"},
    {0x2081, "oneinferior"},
    {0x2082, "twoinferior"},
    {0x2083, "threeinferior"},
    {0x2084, "fourinferior"},
    {0x2085, "fiveinferior"},
    {0x2086, "sixinferior"},
    {0x2087, "seveninferior"},
    {0x2088, "eightinferior"},
    {0x2089, "nineinferior"},
    {0x208D, "parenleftinferior"},
    {0x208E, "parenrightinferior"},
    {0x20A1, "colonmonetary"},
    {0x20A3, "franc"},
    {0x20A4, "lira"},
    {0x20A7, "peseta"},
    {0x20AA, "afii57636"},
    {0x20AB, "dong"},
    {0x20AC, "Euro"},
    {0x2105, "afii61248"},
    {0x2111, "Ifraktur"},
    {0x2113, "afii61289"},
    {0x2116, "afii61352"},
    {0x2118, "weierstrass"},
    {0x211C, "Rfraktur"},
    {0x211E, "prescription"},
    {0x2122, "trademark"},
    {0x2126, "Omega"},
    {0x212E, "estimated"},
    {0x2135, "aleph"},
    {0x2153, "onethird"},
    {0x2154, "twothirds"},
    {0x215B, "oneeighth"},
    {0x215C, "threeeighths"},
    {0x215D, "fiveeighths"},
    {0x215E, "seveneighths"},
    {0x2190, "arrowleft"},
    {0x2191, "arrowup"},
    {0x2192, "arrowright"},
    {0x2193, "arrowdown"},
    {0x2194, "arrowboth"},
    {0x2195, "arrowupdn"},
    {0x21A8, "arrowupdnbse"},
    {0x21B5, "carriagereturn"},
    {0x21D0, "arrowdblleft"},
    {0x21D1, "arrowdblup"},
    {0x21D2, "arrowdblright"},
    {0x21D3, "arrowdbldown"},
    {0x21D4, "arrowdblboth"},
    {0x2200, "universal"},
    {0x2202, "partialdiff"},
    {0x2203, "existential"},
    {0x2205, "emptyset"},
    {0x2206, "Delta"},
    {0x2207, "gradient"},
    {0x2208, "element"},
    {0x2209, "notelement"},
    {0x220B, "suchthat"},
    {0x220F, "product"},
    {0x2211, "summation"},
    {0x2212, "minus"},
    {0x2215, "fraction"},
    {0x2217, "asteriskmath"},
    {0x2219, "periodcentered"},
    {0x221A, "radical"},
    {0x221D, "proportional"},
    {0x221E, "infinity"},
    {0x221F, "orthogonal"},
    {0x2220, "angle"},
    {0x2227, "logicaland"},
    {0x2228, "logicalor"},
    {0x2229, "intersection"},
    {0x222A, "union"},
    {0x222B, "integral"},
    {0x2234, "therefore"},
    {0x223C, "similar"},
    {0x2245, "congruent"},
    {0x2248, "approxequal"},
    {0x2260, "notequal"},
    {0x2261, "equivalence"},
    {0x2264, "lessequal"},
    {0x2265, "greaterequal"},
    {0x2282, "propersubset"},
    {0x2283, "propersuperset"},
    {0x2284, "notsubset"},
    {0x2286, "reflexsubset"},
    {0x2287, "reflexsuperset"},
    {0x2295, "circleplus"},
    {0x2297, "circlemultiply"},
    {0x22A5, "perpendicular"},
    {0x22C5, "dotmath"},
    {0x2302, "house"},
    {0x2310, "revlogicalnot"},
    {0x2320, "integraltp"},
    {0x2321, "integralbt"},
    {0x2329, "angleleft"},
    {0x232A, "angleright"},
    {0x2500, "SF100000"},
    {0x2502, "SF110000"},
    {0x250C, "SF010000"},
    {0x2510, "SF030000"},
    {0x2514, "SF020000"},
    {0x2518, "SF040000"},
    {0x251C, "SF080000"},
    {0x2524, "SF090000"},
    {0x252C, "SF060000"},
    {0x2534, "SF070000"},
    {0x253C, "SF050000"},
    {0x2550, "SF430000"},
    {0x2551, "SF240000"},
    {0x2552, "SF510000"},
    {0x2553, "SF520000"},
    {0x2554, "SF390000"},
    {0x2555, "SF220000"},
    {0x2556, "SF210000"},
    {0x2557, "SF250000"},
    {0x2558, "SF500000"},
    {0x2559, "SF490000"},
    {0x255A, "SF380000"},
    {0x255B, "SF280000"},
    {0x255C, "SF270000"},
    {0x255D, "SF260000"},
    {0x255E, "SF360000"},
    {0x255F, "SF370000"},
    {0x2560, "SF420000"},
    {0x2561, "SF190000"},
    {0x2562, "SF200000"},
    {0x2563, "SF230000"},
    {0x2564, "SF470000"},
    {0x2565, "SF480000"},
    {0x2566, "SF410000"},
    {0x2567, "SF450000"},
    {0x2568, "SF460000"},
    {0x2569, "SF400000"},
    {0x256A, "SF540000"},
    {0x256B, "SF530000"},
    {0x256C, "SF440000"},
    {0x2580, "upblock"},
    {0x2584, "dnblock"},
    {0x2588, "block"},
    {0x258C, "lfblock"},
    {0x2590, "rtblock"},
    {0x2591, "ltshade"},
    {0x2592, "shade"},
    {0x2593, "dkshade"},
    {0x25A0, "filledbox"},
    {0x25A1, "H22073"},
    {0x25AA, "H18543"},
    {0x25AB, "H18551"},
    {0x25AC, "filledrect"},
    {0x25B2, "triagup"},
    {0x25BA, "triagrt"},
    {0x25BC, "triagdn"},
    {0x25C4, "triaglf"},
    {0x25CA, "lozenge"},
    {0x25CB, "circle"},
    {0x25CF, "H18533"},
    {0x25D8, "invbullet"},
    {0x25D9, "invcircle"},
    {0x25E6, "openbullet"},
    {0x263A, "smileface"},
    {0x263B, "invsmileface"},
    {0x263C, "sun"},
    {0x2640, "female"},
    {0x2642, "male"},
    {0x2660, "spade"},
    {0x2663, "club"},
    {0x2665, "heart"},
    {0x2666, "diamond"},
    {0x266A, "musicalnote"},
    {0x266B, "musicalnotedbl"},
    {0xF6BE, "dotlessj"},
    {0xF6BF, "LL"},
    {0xF6C0, "ll"},
    {0xF6C1, "Scedilla"},
    {0xF6C2, "scedilla"},
    {0xF6C3, "commaaccent"},
    {0xF6C4, "afii10063"},
    {0xF6C5, "afii10064"},
    {0xF6C6, "afii10192"},
    {0xF6C7, "afii10831"},
    {0xF6C8, "afii10832"},
    {0xF6C9, "Acute"},
    {0xF6CA, "Caron"},
    {0xF6CB, "Dieresis"},
    {0xF6CC, "DieresisAcute"},
    {0xF6CD, "DieresisGrave"},
    {0xF6CE, "Grave"},
    {0xF6CF, "Hungarumlaut"},
    {0xF6D0, "Macron"},
    {0xF6D1, "cyrBreve"},
    {0xF6D2, "cyrFlex"},
    {0xF6D3, "dblGrave"},
    {0xF6D4, "cyrbreve"},
    {0xF6D5, "cyrflex"},
    {0xF6D6, "dblgrave"},
    {0xF6D7, "dieresisacute"},
    {0xF6D8, "dieresisgrave"},
    {0xF6D9, "copyrightserif"},
    {0xF6DA, "registerserif"},
    {0xF6DB, "trademarkserif"},
    {0xF6DC, "onefitted"},
    {0xF6DD, "rupiah"},
    {0xF6DE, "threequartersemdash"},
    {0xF6DF, "centinferior"},
    {0xF6E0, "centsuperior"},
    {0xF6E1, "commainferior"},
    {0xF6E2, "commasuperior"},
    {0xF6E3, "dollarinferior"},
    {0xF6E4, "dollarsuperior"},
    {0xF6E5, "hypheninferior"},
    {0xF6E6, "hyphensuperior"},
    {0xF6E7, "periodinferior"},
    {0xF6E8, "periodsuperior"},
    {0xF6E9, "asuperior"},
    {0xF6EA, "bsuperior"},
    {0xF6EB, "dsuperior"},
    {0xF6EC, "esuperior"},
    {0xF6ED, "isuperior"},
    {0xF6EE, "lsuperior"},
    {0xF6EF, "msuperior"},
    {0xF6F0, "osuperior"},
    {0xF6F1, "rsuperior"},
    {0xF6F2, "ssuperior"},
    {0xF6F3, "tsuperior"},
    {0xF6F4, "Brevesmall"},
    {0xF6F5, "Caronsmall"},
    {0xF6F6, "Circumflexsmall"},
    {0xF6F7, "Dotaccentsmall"},
    {0xF6F8, "Hungarumlautsmall"},
    {0xF6F9, "Lslashsmall"},
    {0xF6FA, "OEsmall"},
    {0xF6FB, "Ogoneksmall"},
    {0xF6FC, "Ringsmall"},
    {0xF6FD, "Scaronsmall"},
    {0xF6FE, "Tildesmall"},
    {0xF6FF, "Zcaronsmall"},
    {0xF721, "exclamsmall"},
    {0xF724, "dollaroldstyle"},
    {0xF726, "ampersandsmall"},
    {0xF730, "zerooldstyle"},
    {0xF731, "oneoldstyle"},
    {0xF732, "twooldstyle"},
    {0xF733, "threeoldstyle"},
    {0xF734, "fouroldstyle"},
    {0xF735, "fiveoldstyle"},
    {0xF736, "sixoldstyle"},
    {0xF737, "sevenoldstyle"},
    {0xF738, "eightoldstyle"},
    {0xF739, "nineoldstyle"},
    {0xF73F, "questionsmall"},
    {0xF760, "Gravesmall"},
    {0xF761, "Asmall"},
    {0xF762, "Bsmall"},
    {0xF763, "Csmall"},
    {0xF764, "Dsmall"},
    {0xF765, "Esmall"},
    {0xF766, "Fsmall"},
    {0xF767, "Gsmall"},
    {0xF768, "Hsmall"},
    {0xF769, "Ismall"},
    {0xF76A, "Jsmall"},
    {0xF76B, "Ksmall"},
    {0xF76C, "Lsmall"},
    {0xF76D, "Msmall"},
    {0xF76E, "Nsmall"},
    {0xF76F, "Osmall"},
    {0xF770, "Psmall"},
    {0xF771, "Qsmall"},
    {0xF772, "Rsmall"},
    {0xF773, "Ssmall"},
    {0xF774, "Tsmall"},
    {0xF775, "Usmall"},
    {0xF776, "Vsmall"},
    {0xF777, "Wsmall"},
    {0xF778, "Xsmall"},
    {0xF779, "Ysmall"},
    {0xF77A, "Zsmall"},
    {0xF7A1, "exclamdownsmall"},
    {0xF7A2, "centoldstyle"},
    {0xF7A8, "Dieresissmall"},
    {0xF7AF, "Macronsmall"},
    {0xF7B4, "Acutesmall"},
    {0xF7B8, "Cedillasmall"},
    {0xF7BF, "questiondownsmall"},
    {0xF7E0, "Agravesmall"},
    {0xF7E1, "Aacutesmall"},
    {0xF7E2, "Acircumflexsmall"},
    {0xF7E3, "Atildesmall"},
    {0xF7E4, "Adieresissmall"},
    {0xF7E5, "Aringsmall"},
    {0xF7E6, "AEsmall"},
    {0xF7E7, "Ccedillasmall"},
    {0xF7E8, "Egravesmall"},
    {0xF7E9, "Eacutesmall"},
    {0xF7EA, "Ecircumflexsmall"},
    {0xF7EB, "Edieresissmall"},
    {0xF7EC, "Igravesmall"},
    {0xF7ED, "Iacutesmall"},
    {0xF7EE, "Icircumflexsmall"},
    {0xF7EF, "Idieresissmall"},
    {0xF7F0, "Ethsmall"},
    {0xF7F1, "Ntildesmall"},
    {0xF7F2, "Ogravesmall"},
    {0xF7F3, "Oacutesmall"},
    {0xF7F4, "Ocircumflexsmall"},
    {0xF7F5, "Otildesmall"},
    {0xF7F6, "Odieresissmall"},
    {0xF7F8, "Oslashsmall"},
    {0xF7F9, "Ugravesmall"},
    {0xF7FA, "Uacutesmall"},
    {0xF7FB, "Ucircumflexsmall"},
    {0xF7FC, "Udieresissmall"},
    {0xF7FD, "Yacutesmall"},
    {0xF7FE, "Thornsmall"},
    {0xF7FF, "Ydieresissmall"},
    {0xF8E5, "radicalex"},
    {0xF8E6, "arrowvertex"},
    {0xF8E7, "arrowhorizex"},
    {0xF8E8, "registersans"},
    {0xF8E9, "copyrightsans"},
    {0xF8EA, "trademarksans"},
    {0xF8EB, "parenlefttp"},
    {0xF8EC, "parenleftex"},
    {0xF8ED, "parenleftbt"},
    {0xF8EE, "bracketlefttp"},
    {0xF8EF, "bracketleftex"},
    {0xF8F0, "bracketleftbt"},
    {0xF8F1, "bracelefttp"},
    {0xF8F2, "braceleftmid"},
    {0xF8F3, "braceleftbt"},
    {0xF8F4, "braceex"},
    {0xF8F5, "integralex"},
    {0xF8F6, "parenrighttp"},
    {0xF8F7, "parenrightex"},
    {0xF8F8, "parenrightbt"},
    {0xF8F9, "bracketrighttp"},
    {0xF8FA, "bracketrightex"},
    {0xF8FB, "bracketrightbt"},
    {0xF8FC, "bracerighttp"},
    {0xF8FD, "bracerightmid"},
    {0xF8FE, "bracerightbt"},
    {0xFB00, "ff"},
    {0xFB01, "fi"},
    {0xFB02, "fl"},
    {0xFB03, "ffi"},
    {0xFB04, "ffl"},
    {0xFB1F, "afii57705"},
    {0xFB2A, "afii57694"},
    {0xFB2B, "afii57695"},
    {0xFB35, "afii57723"},
    {0xFB4B, "afii57700"},
    {0xFFFF, NULL}
};

const char*
PdfEncodingDef::UnicodeToGryphName(unsigned short unicode) 
{
    static const char* NOTDEF = PDF_CHAR_NOTDEF;
    const pdf_unicode_gryph_pair* map = UNICODE_GRYPH_NAME_MAP;

    while (map->unicode <= unicode) {
        if (map->unicode == unicode)
            return map->gryph_name;
        map++;
    }
    return NOTDEF;
}

unsigned short
PdfEncodingDef::GryphNameToUnicode(const char* gryph_name)
{
    const pdf_unicode_gryph_pair* map = UNICODE_GRYPH_NAME_MAP;

    while (map->unicode != 0xFFFF) {
        if (strcmp(gryph_name, map->gryph_name) == 0)
            return map->unicode;
        map++;
    }
    return 0x0000;
}

const char*
PdfEncodingDef::PdfEncodingToString(pdf_encoding encoding)
{
    int idx = -1;
    for (int i = 0; i < (int)PDF_ENCODING_EOF; i++){
        if (pdf_encoding_names[i].encode == encoding) {
            idx = i;
            break;
        }
    }

    if (idx < 0) {
        PDF_DEBUG_PRINT(("ERROR: internal error invalid encoding [%d].\n",
             encoding));
        return NULL;
    }

    return pdf_encoding_names[idx].encode_name;
}

pdf_encoding
PdfEncodingDef::StringToPdfEncoding(const char* encoding)
{
    if (encoding == NULL)
        return PDF_ENCODING_EOF;

    for (int i = 0; i < (int)PDF_ENCODING_EOF; i++){
        if (strcmp(pdf_encoding_names[i].encode_name, encoding) == 0)
            return (pdf_encoding)i;
    }

    if (strcmp(encoding, "AdobeStandardEncoding") == 0)
        return PDF_STANDARD_ENCODING;
    return PDF_ENCODING_EOF;
}

/*----------------------------------------------------------------------------*/
/*------ PdfStandardEncoding -------------------------------------------------*/

PdfStandardEncoding::PdfStandardEncoding()
        : PdfEncodingDef()
{
    static const unsigned int FIRST_CHAR= 32;
    static const unsigned int LAST_CHAR= 255;

    const unsigned short UNICODE_MAP_STANDARD[] = {
        0x0020, 0x0021, 0x0022, 0x0023, 0x0024, 0x0025, 0x0026, 0x2019,
        0x0028, 0x0029, 0x002A, 0x002B, 0x002C, 0x002D, 0x002E, 0x002F,
        0x0030, 0x0031, 0x0032, 0x0033, 0x0034, 0x0035, 0x0036, 0x0037,
        0x0038, 0x0039, 0x003A, 0x003B, 0x003C, 0x003D, 0x003E, 0x003F,
        0x0040, 0x0041, 0x0042, 0x0043, 0x0044, 0x0045, 0x0046, 0x0047,
        0x0048, 0x0049, 0x004A, 0x004B, 0x004C, 0x004D, 0x004E, 0x004F,
        0x0050, 0x0051, 0x0052, 0x0053, 0x0054, 0x0055, 0x0056, 0x0057,
        0x0058, 0x0059, 0x005A, 0x005B, 0x005C, 0x005D, 0x005E, 0x005F,
        0x2018, 0x0061, 0x0062, 0x0063, 0x0064, 0x0065, 0x0066, 0x0067,
        0x0068, 0x0069, 0x006A, 0x006B, 0x006C, 0x006D, 0x006E, 0x006F,
        0x0070, 0x0071, 0x0072, 0x0073, 0x0074, 0x0075, 0x0076, 0x0077,
        0x0078, 0x0079, 0x007A, 0x007B, 0x007C, 0x007D, 0x007E, 0x0000,
        0x0000, 0x0000, 0x0000, 0x0000, 0x00D1, 0x0000, 0x0000, 0x0000,
        0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0000, 0x00A1, 0x00A2, 0x00A3, 0x2044, 0x00A5, 0x0192, 0x00A7,
        0x00A4, 0x0027, 0x201C, 0x00AB, 0x2039, 0x203A, 0xFB01, 0xFB02,
        0x0000, 0x2013, 0x2020, 0x2021, 0x00B7, 0x0000, 0x00B6, 0x2022,
        0x201A, 0x201E, 0x201D, 0x00BB, 0x2026, 0x2030, 0x0000, 0x00BF,
        0x0000, 0x0060, 0x00B4, 0x02C6, 0x02DC, 0x00AF, 0x02D8, 0x02D9,
        0x00A8, 0x0000, 0x02DA, 0x00B8, 0x0000, 0x02DD, 0x02DB, 0x02C7,
        0x2014, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0000, 0x00C6, 0x0000, 0x00AA, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0141, 0x00D8, 0x0152, 0x00BA, 0x0000, 0x0000, 0x0000, 0x0000,
        0x0000, 0x00E6, 0x0000, 0x0000, 0x0000, 0x0131, 0x0000, 0x0000,
        0x0142, 0x00F8, 0x0153, 0x00DF, 0x0000, 0x0000, 0x0000, 0x0000
    };

    SetUnicodeArray(FIRST_CHAR, LAST_CHAR, UNICODE_MAP_STANDARD);
    fBaseEncoding = PDF_STANDARD_ENCODING;
}

/*----------------------------------------------------------------------------*/
/*------ WinAnsiEncoding -----------------------------------------------------*/

PdfWinAnsiEncoding::PdfWinAnsiEncoding()
        : PdfEncodingDef()
{
    static const unsigned int FIRST_CHAR= 32;
    static const unsigned int LAST_CHAR= 255;
    
    static const unsigned short UNICODE_MAP_WIN_ANSI[224] = {
        0x0020, 0x0021, 0x0022, 0x0023, 0x0024, 0x0025, 0x0026, 0x0027, 
        0x0028, 0x0029, 0x002A, 0x002B, 0x002C, 0x002D, 0x002E, 0x002F, 
        0x0030, 0x0031, 0x0032, 0x0033, 0x0034, 0x0035, 0x0036, 0x0037, 
        0x0038, 0x0039, 0x003A, 0x003B, 0x003C, 0x003D, 0x003E, 0x003F, 
        0x0040, 0x0041, 0x0042, 0x0043, 0x0044, 0x0045, 0x0046, 0x0047, 
        0x0048, 0x0049, 0x004A, 0x004B, 0x004C, 0x004D, 0x004E, 0x004F, 
        0x0050, 0x0051, 0x0052, 0x0053, 0x0054, 0x0055, 0x0056, 0x0057, 
        0x0058, 0x0059, 0x005A, 0x005B, 0x005C, 0x005D, 0x005E, 0x005F, 
        0x0060, 0x0061, 0x0062, 0x0063, 0x0064, 0x0065, 0x0066, 0x0067, 
        0x0068, 0x0069, 0x006A, 0x006B, 0x006C, 0x006D, 0x006E, 0x006F, 
        0x0070, 0x0071, 0x0072, 0x0073, 0x0074, 0x0075, 0x0076, 0x0077, 
        0x0078, 0x0079, 0x007A, 0x007B, 0x007C, 0x007D, 0x007E, 0x0000, 
        0x20AC, 0x0000, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 
        0x02C6, 0x2030, 0x0160, 0x2039, 0x0152, 0x0000, 0x017D, 0x0000, 
        0x0000, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022, 0x2013, 0x2014, 
        0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x0000, 0x017E, 0x0178, 
        0x0000, 0x00A1, 0x00A2, 0x00A3, 0x00A4, 0x00A5, 0x00A6, 0x00A7, 
        0x00A8, 0x00A9, 0x00AA, 0x00AB, 0x00AC, 0x0000, 0x00AE, 0x00AF, 
        0x02DA, 0x00B1, 0x00B2, 0x00B3, 0x00B4, 0x00B5, 0x00B6, 0x00B7, 
        0x00B8, 0x00B9, 0x00BA, 0x00BB, 0x00BC, 0x00BD, 0x00BE, 0x00BF, 
        0x00C0, 0x00C1, 0x00C2, 0x00C3, 0x00C4, 0x00C5, 0x00C6, 0x00C7, 
        0x00C8, 0x00C9, 0x00CA, 0x00CB, 0x00CC, 0x00CD, 0x00CE, 0x00CF, 
        0x00D0, 0x00D1, 0x00D2, 0x00D3, 0x00D4, 0x00D5, 0x00D6, 0x00D7, 
        0x00D8, 0x00D9, 0x00DA, 0x00DB, 0x00DC, 0x00DD, 0x00DE, 0x00DF, 
        0x00E0, 0x00E1, 0x00E2, 0x00E3, 0x00E4, 0x00E5, 0x00E6, 0x00E7, 
        0x00E8, 0x00E9, 0x00EA, 0x00EB, 0x00EC, 0x00ED, 0x00EE, 0x00EF, 
        0x00F0, 0x00F1, 0x00F2, 0x00F3, 0x00F4, 0x00F5, 0x00F6, 0x00F7, 
        0x00F8, 0x00F9, 0x00FA, 0x00FB, 0x00FC, 0x00FD, 0x00FE, 0x00FF
    };

    SetUnicodeArray(FIRST_CHAR, LAST_CHAR, UNICODE_MAP_WIN_ANSI);
    fBaseEncoding = PDF_WIN_ANSI_ENCODING;
}

/*----------------------------------------------------------------------------*/
/*------ PdfMacRomanEncoding -------------------------------------------------*/

PdfMacRomanEncoding::PdfMacRomanEncoding()
        : PdfEncodingDef()
{
    static const unsigned int FIRST_CHAR= 32;
    static const unsigned int LAST_CHAR= 255;
    
    static const unsigned short UNICODE_MAP_MAC_ROMAN[224] = {
        0x0020, 0x0021, 0x0022, 0x0023, 0x0024, 0x0025, 0x0026, 0x0027,
        0x0028, 0x0029, 0x002A, 0x002B, 0x002C, 0x002D, 0x002E, 0x002F,
        0x0030, 0x0031, 0x0032, 0x0033, 0x0034, 0x0035, 0x0036, 0x0037,
        0x0038, 0x0039, 0x003A, 0x003B, 0x003C, 0x003D, 0x003E, 0x003F,
        0x0040, 0x0041, 0x0042, 0x0043, 0x0044, 0x0045, 0x0046, 0x0047,
        0x0048, 0x0049, 0x004A, 0x004B, 0x004C, 0x004D, 0x004E, 0x004F,
        0x0050, 0x0051, 0x0052, 0x0053, 0x0054, 0x0055, 0x0056, 0x0057,
        0x0058, 0x0059, 0x005A, 0x005B, 0x005C, 0x005D, 0x005E, 0x005F,
        0x0060, 0x0061, 0x0062, 0x0063, 0x0064, 0x0065, 0x0066, 0x0067,
        0x0068, 0x0069, 0x006A, 0x006B, 0x006C, 0x006D, 0x006E, 0x006F,
        0x0070, 0x0071, 0x0072, 0x0073, 0x0074, 0x0075, 0x0076, 0x0077,
        0x0078, 0x0079, 0x007A, 0x007B, 0x007C, 0x007D, 0x007E, 0x0000,
        0x00C4, 0x00C5, 0x00C7, 0x00C9, 0x0000, 0x00D6, 0x00DC, 0x00E1,
        0x00E0, 0x00E2, 0x00E4, 0x00E3, 0x00E5, 0x00E7, 0x00E9, 0x00E8,
        0x00EA, 0x00EB, 0x00ED, 0x00EC, 0x00EE, 0x00EF, 0x00F1, 0x00F3,
        0x00F2, 0x00F4, 0x00F6, 0x00F5, 0x00FA, 0x00F9, 0x00FB, 0x00FC,
        0x2020, 0x00B0, 0x00A2, 0x00A3, 0x00A7, 0x2022, 0x00B6, 0x00DF,
        0x00AE, 0x00A9, 0x2122, 0x00B4, 0x00A8, 0x0000, 0x00C6, 0x00D8,
        0x0000, 0x00B1, 0x0000, 0x0000, 0x00A5, 0x00B5, 0x0000, 0x0000,
        0x0000, 0x0000, 0x0000, 0x00AA, 0x00BA, 0x0000, 0x00E6, 0x00F8,
        0x00BF, 0x00A1, 0x00AC, 0x0000, 0x0192, 0x0000, 0x0000, 0x00AB,
        0x00BB, 0x2026, 0x0020, 0x00C0, 0x00C3, 0x00D5, 0x0152, 0x0153,
        0x2013, 0x2014, 0x201C, 0x201D, 0x2018, 0x2019, 0x00F7, 0x0000,
        0x00FF, 0x0178, 0x2044, 0x00A4, 0x2039, 0x203A, 0xFB01, 0xFB02,
        0x2021, 0x00B7, 0x201A, 0x201E, 0x2030, 0x00C2, 0x00CA, 0x00C1,
        0x00CB, 0x00C8, 0x00CD, 0x00CE, 0x00CF, 0x00CC, 0x00D3, 0x00D4,
        0x0000, 0x00D2, 0x00DA, 0x00DB, 0x00D9, 0x0131, 0x02C6, 0x02DC,
        0x00AF, 0x02D8, 0x02D9, 0x02DA, 0x00B8, 0x02DD, 0x02DB, 0x02C7
    };

    SetUnicodeArray(FIRST_CHAR, LAST_CHAR, UNICODE_MAP_MAC_ROMAN);
    fBaseEncoding = PDF_MAC_ROMAN_ENCODING;
}

/*----------------------------------------------------------------------------*/

