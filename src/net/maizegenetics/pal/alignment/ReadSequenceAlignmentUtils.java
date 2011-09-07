// ReadSequenceAlignmentUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.io.FormattedInput;
import net.maizegenetics.pal.io.InputSource;

import java.io.IOException;
import java.io.PushbackReader;
import java.util.Vector;
import net.maizegenetics.pal.datatype.IUPACNucleotides;

/**
 * reads aligned sequence data from plain text files.<p>
 *
 * recognizes PHYLIP 3.4 INTERLEAVED,
 *              PHYLIP SEQUENTIAL,
 *              CLUSTAL and derived formats.<p>
 *
 * Other features:
 * - the dot as "copy character" is recognized,
 * - all base characters are capitalized,
 * - automatic data type estimation
 * - determination of corresponding base frequencies.
 *
 * @version $Id: ReadSequenceAlignmentUtils.java,v 1.1 2009/07/23 19:39:41 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class ReadSequenceAlignmentUtils {

    private ReadSequenceAlignmentUtils() {
        // to prevent instantiation of this utility class.
    }

    /** read from stream, while letting the user set the maximum label (taxa) label length
     * as non-standard phylip formats now exist*/
    public static Alignment readBasicAlignments(PushbackReader input, int maxLabelLength)
            throws IOException {
        Alignment saa = null;
        saa = readFile(input, maxLabelLength);
        return saa;
    }

    /** read from file, while letting the user set the maximum label (taxa) label length
     * as non-standard phylip formats now exist*/
    public static Alignment readBasicAlignments(String file, int maxLabelLength)
            throws IOException {
        PushbackReader input = InputSource.openFile(file);
        Alignment saa = readBasicAlignments(input, maxLabelLength);
        input.close();
        return saa;
    }

    private static final boolean isType(PushbackReader in, String id) throws IOException {
        FormattedInput fi = FormattedInput.getInstance();
        for (int i = 0; i < id.length(); i++) {
            int c = fi.readNextChar(in);
            if (c != id.charAt(i)) {
                in.unread(c);
                return false;
            }
        }
        return true;
    }

    private static Alignment readFile(PushbackReader in, int maxLabelLength) throws IOException {
        FormattedInput fi = FormattedInput.getInstance();
        Alignment saa = null;
        if (isType(in, "CLUSTAL")) {
            fi.nextLine(in);
            saa = readCLUSTALW(in, maxLabelLength);
        } else if (isType(in, "#NEXUS") || isType(in, "#nexus")) {
            saa = readNEXUS(in);
        } else {
            saa = readPHYLIP(in, maxLabelLength);
        }
        return saa;
    }

    private static Alignment readNEXUS(PushbackReader in) {
        FormattedInput fi = FormattedInput.getInstance();
        IdGroup idGroup;
        int numSeqs = 0, numSites = 0;
        char[][] data = null;
        int c, seq = 0, pos = 0;
        try {
            // Find start block
            // going to try to find the beggining of the data
            // Let's find information about the file in a massive way. This is the
            // best way I've thought of so far. :)
            boolean seqsOK = false;
            boolean sitesOK = false;
            String tempWord = null;

            tempWord = fi.readWord(in);
            while (!(seqsOK && sitesOK)) {
                if (tempWord.length() >= 5) {
                    if ((tempWord.substring(0, 5).equals("nchar") || tempWord.substring(0, 5).equals("NCHAR"))) {
                        int stringSize = tempWord.length();
                        if (tempWord.substring(stringSize - 1).equals(";")) {
                            stringSize = stringSize - 1;
                        }
                        numSites = new Integer(tempWord.substring(6, stringSize)).intValue();
                        sitesOK = true;
                    } else if ((tempWord.substring(0, 4).equals("ntax") || tempWord.substring(0, 4).equals("NTAX"))) {
                        int stringSize = tempWord.length();
                        if (tempWord.substring(stringSize - 1).equals(";")) {
                            stringSize = stringSize - 1;
                        }
                        numSeqs = new Integer(tempWord.substring(5, stringSize)).intValue();
                        seqsOK = true;
                    }
                }
                //String ahpois = tempWord.substring(0,5) ;
                tempWord = fi.readWord(in);
            }


            // Now let's look for the initial of the matrix :-) We first have to find the "format" line
            //tempWord = fi.readWord(in) ;
            while (!(tempWord.equals("Format") || tempWord.equals("format") || tempWord.equals("FORMAT"))) {
                fi.nextLine(in);
                tempWord = fi.readWord(in);
            }

            tempWord = fi.readWord(in);
            while (!(tempWord.equals("Matrix") || tempWord.equals("matrix") || tempWord.equals("MATRIX"))) {
                fi.nextLine(in);
                tempWord = fi.readWord(in);
            }

            //Reserve memory (as per readPHYLIP)
            String[] identifiers = new String[numSeqs];
            data = new char[numSeqs][numSites];

            // And now back to our regular show.. :-)
            //c = fi.readNextChar(in);
            //in.unread(c);
            for (seq = 0; seq < numSeqs; seq++) {
                // Go to next block
                c = fi.readNextChar(in);
                in.unread(c);

                // Read label
                identifiers[seq] = fi.readLabel(in, 30).toUpperCase();

                // Read sequences
                for (pos = 0; pos < numSites; pos++) {
                    data[seq][pos] = (char) fi.readNextChar(in);

                    if (data[0][pos] == '.') {
                        if (seq == 0) {
                            throw new IllegalArgumentException(
                                    "Copy character (.) in first sequence not allowed (pos. " + (pos + 1) + ")");
                        } else {
                            data[seq][pos] = data[0][pos];
                        }
                    }
                }
            }

            idGroup = new SimpleIdGroup(identifiers);
        } catch (IOException e) {
            throw new IllegalArgumentException("IO error after pos. " + (pos + 1) + ", seq. " + (seq + 1));
        }
        String[] s = new String[numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            s[i] = (new String(data[i])).toUpperCase();
        }
        //SimpleAlignment saa = new SimpleAlignment(idGroup, s, new Nucleotides());
//        SimpleAlignment saa = SimpleAlignment.getInstance(idGroup, s, new Nucleotides());
        String[] sites = new String[numSites];
        int[] positions = new int[numSites];
        for (int i = 0; i < numSites; i++) {
        	positions[i] = i;
        	sites[i] = Integer.toString(i);
        }
        SimpleAlignment saa = new SimpleAlignment(idGroup, s, new IUPACNucleotides(), null, positions, new byte[0], "", new float[0][0], sites);

        return saa;
    }

    // Read alignment (in CLUSTAL W format)
    private static Alignment readCLUSTALW(PushbackReader in, int maxLabelLength) {
        FormattedInput fi = FormattedInput.getInstance();
        IdGroup idGroup;
        int numSeqs = 0, numSites = 0, lineLength = 0;
        char[][] data = null;
        Vector names, seqs, sites;
        int EOF = -1;
        int c, seq = 0, pos = 0;
        try {
            // Find start block
            c = fi.readNextChar(in);
            in.unread(c);
            names = new Vector();
            seqs = new Vector();
            sites = new Vector();

            // Reading first data block
            c = in.read();
            while (!Character.isWhitespace((char) c)) {
                in.unread(c);
                names.addElement(fi.readLabel(in, maxLabelLength));
                lineLength = readSeqLineC(in, seq, pos, sites, seqs, fi, lineLength);
                seq++;
                c = in.read();
            }
            in.unread(c);
            // Skip CLUSTAL W status line
            fi.nextLine(in);

            pos += lineLength;
            numSeqs = seq;

            // Go to next block
            c = fi.readNextChar(in);
            in.unread(c);

            // Reading remaining blocks
            while (c != EOF) {
                for (seq = 0; seq < numSeqs; seq++) {
                    // goto next blank
                    do {
                        c = in.read();
                        if (c < 0) {
                            throw new IllegalArgumentException("Unexpected end of file exception!");
                        }
                    } while (Character.isWhitespace((char) c));
                    lineLength = readSeqLineC(in, seq, pos, sites, seqs, fi, lineLength);

                }

                // Skip CLUSTAL W status line
                fi.nextLine(in);

                pos += lineLength;

                // Go to next block
                c = fi.readNextChar(in);
                in.unread(c);
            }

            numSites = pos;

            // Copy to array
            String[] identifiers = new String[numSeqs];
            data = new char[numSeqs][numSites];
            for (int i = 0; i < numSeqs; i++) {
                identifiers[i] = (String) names.elementAt(i);
            }
            idGroup = new SimpleIdGroup(identifiers);
            for (int i = 0; i < numSeqs; i++) {
                for (int j = 0; j < numSites; j++) {
                    data[i][j] =
                            ((Character) ((Vector) seqs.elementAt(i)).elementAt(j)).charValue();
                }
            }

            // Help garbage collector
            names = null;
            for (int i = 0; i < numSeqs; i++) {
                ((Vector) seqs.elementAt(i)).removeAllElements();
            }
            seqs = null;
        } catch (IOException e) {
            throw new IllegalArgumentException("IO error after pos. " + (pos + 1) + ", seq. " + (seq + 1));
        }
        String[] s = new String[numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            s[i] = (new String(data[i])).toUpperCase();
        }
//        SimpleAlignment saa = SimpleAlignment.getInstance(idGroup, s, new Nucleotides());
        String[] siteNames = new String[numSites];
        int[] positions = new int[numSites];
        for (int i = 0; i < numSites; i++) {
        	positions[i] = i;
        	siteNames[i] = Integer.toString(i);
        }
        SimpleAlignment saa = new SimpleAlignment(idGroup, s, new IUPACNucleotides(), null, positions, new byte[0], "", new float[0][0], siteNames);

        return saa;
    }

    private static int readSeqLineC(PushbackReader in, int s, int pos, Vector sites, Vector seqs, FormattedInput fi, int lineLength)
            throws IOException {
        int c;
        if (pos == 0) {
            sites = new Vector();
            seqs.addElement(sites);
        } else {
            sites = (Vector) seqs.elementAt(s);
        }

        if (s == 0) {
            String thisLine = fi.readLine(in, false);
            lineLength = thisLine.length();

            for (int i = 0; i < lineLength; i++) {
                c = thisLine.charAt(i);
                if (c == '.') {
                    throw new IllegalArgumentException("Copy character (.) in first sequence not allowed (pos. " + (i + pos + 1) + ")");
                }
                sites.addElement(new Character((char) c));
            }
        } else {
            for (int i = 0; i < lineLength; i++) {
                c = fi.readNextChar(in);
                if (c == '.') {
                    c = ((Character) ((Vector) seqs.elementAt(0)).elementAt(pos + i)).charValue();
                }
                sites.addElement(new Character((char) c));
            }
            fi.nextLine(in);
        }
        return lineLength;
    }

    // Read alignment (in PHYLIP 3.4 INTERLEAVED or PHYLIP SEQUENTIAL format)
    private static Alignment readPHYLIP(PushbackReader in, int maxLabelLength) {
        FormattedInput fi = FormattedInput.getInstance();
        IdGroup idGroup;
        int numSeqs = 0, numSites = 0, lineLength = 0;
        char[][] data = null;
        int c, pos = 0, seq = 0;

        try {
            // Parse PHYLIP header line
            numSeqs = fi.readInt(in);
            numSites = fi.readInt(in);

            String[] identifiers = new String[numSeqs];
            data = new char[numSeqs][numSites];


            // Determine whether sequences are in INTERLEAVED
            // or in sequential format
            String header = fi.readLine(in, false);

            boolean interleaved = true;

            if (header.length() > 0) {
                if (header.charAt(0) == 'S') {
                    interleaved = false;
                }
            }

            if (interleaved) // PHYLIP INTERLEAVED
            {
    
                // Reading data
                while (pos < numSites) {
                    // Go to next block
                    c = fi.readNextChar(in);
                    in.unread(c);

                    for (seq = 0; seq < numSeqs; seq++) {
                        lineLength = readSeqLineP(in, seq, pos, numSites, data, identifiers, fi, maxLabelLength, lineLength);
                    }
                    pos += lineLength;
                }
            } else // PHYLIP SEQUENTIAL
            {
                //System.out.println("PHYLIP SEQUENTIAL");

                for (seq = 0; seq < numSeqs; seq++) {
                    // Go to next block
                    c = fi.readNextChar(in);
                    in.unread(c);

                    // Read label
                    identifiers[seq] = fi.readLabel(in, maxLabelLength).toUpperCase();

                    // Read sequences
                    for (pos = 0; pos < numSites; pos++) {
                        data[seq][pos] = (char) fi.readNextChar(in);

                        if (data[0][pos] == '.') {
                            if (seq == 0) {
                                throw new IllegalArgumentException(
                                        "Copy character (.) in first sequence not allowed (pos. " + (pos + 1) + ")");
                            } else {
                                data[seq][pos] = data[0][pos];
                            }
                        }
                    }
                }
            }
            idGroup = new SimpleIdGroup(identifiers);
        } catch (IOException e) {
            throw new IllegalArgumentException("IO error after pos. " + (pos + 1) + ", seq. " + (seq + 1));
        }
        String[] s = new String[numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            s[i] = (new String(data[i])).toUpperCase();
        }
//        SimpleAlignment saa = SimpleAlignment.getInstance(idGroup, s, new Nucleotides());
        String[] sites = new String[numSites];
        int[] positions = new int[numSites];
        for (int i = 0; i < numSites; i++) {
        	positions[i] = i;
        	sites[i] = Integer.toString(i);
        }
        SimpleAlignment saa = new SimpleAlignment(idGroup, s, new IUPACNucleotides(), null, positions, new byte[0], "", new float[0][0], sites);
        return saa;
    }

    private static int readSeqLineP(PushbackReader in, int s, int pos, int maxPos, char[][] data, String[] identifiers,
            FormattedInput fi, int maxLabelLength, int lineLength)
            throws IOException {
        if (pos == 0) {
            identifiers[s] = fi.readLabel(in, maxLabelLength).toUpperCase();
        }

        if (s == 0) {
            String thisLine = fi.readLine(in, false);

            if (thisLine.length() > maxPos - pos) {
                lineLength = maxPos - pos;
            } else {
                lineLength = thisLine.length();
            }

            for (int i = 0; i < lineLength; i++) {
                data[0][pos + i] = thisLine.charAt(i);
                if (data[0][pos + i] == '.') {
                    throw new IllegalArgumentException("Copy character (.) in first sequence not allowed (pos. " + (i + pos + 1) + ")");
                }
            }
        } else {
            for (int i = 0; i < lineLength; i++) {
                data[s][pos + i] = (char) fi.readNextChar(in);
                if (data[s][pos + i] == '.') {
                    data[s][pos + i] = data[0][pos + i];
                }
            }
            fi.nextLine(in);
        }
        return lineLength;
    }
}
