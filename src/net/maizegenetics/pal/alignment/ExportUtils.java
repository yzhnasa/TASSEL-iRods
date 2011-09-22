/*
 * ExportUtils
 */
package net.maizegenetics.pal.alignment;

import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintWriter;

import net.maizegenetics.util.Utils;

import java.util.regex.Pattern;

import net.maizegenetics.pal.io.FormattedOutput;

import org.apache.log4j.Logger;

/**
 * The class exports PAL alignment datatypes to
 * various file formats.
 *
 * @author Jon
 */
public class ExportUtils {

    private static final Logger myLogger = Logger.getLogger(ExportUtils.class);
    static FormattedOutput format = FormattedOutput.getInstance();

    private ExportUtils() {
        // Utility Class - do not instantiate.
    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error checking
     * @param alignment array of alignemnts
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static void writeToHapmap(Alignment alignment, boolean diploid, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".hmp.txt")), 1000000);
            bw.write("rs#");
            bw.write(delimChar);
            bw.write("alleles");
            bw.write(delimChar);
            bw.write("chrom");
            bw.write(delimChar);
            bw.write("pos");
            bw.write(delimChar);
            bw.write("strand");
            bw.write(delimChar);
            bw.write("assembly#");
            bw.write(delimChar);
            bw.write("center");
            bw.write(delimChar);
            bw.write("protLSID");
            bw.write(delimChar);
            bw.write("assayLSID");
            bw.write(delimChar);
            bw.write("panelLSID");
            bw.write(delimChar);
            bw.write("QCcode");
            bw.write(delimChar);
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                //not completely sure this does what I want, I need to access the
                //accession name from every alleleBLOB in bytes [52-201] but there
                //doesn't seem to be a method to access that in Alignment
                String sequenceID = alignment.getIdGroup().getIdentifier(taxa).getFullName().trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
                StringBuilder sb = new StringBuilder(alignment.getSNPID(site));
                sb.append(delimChar);
                byte[] alleles = alignment.getAlleles(site);
                int numAlleles = alleles.length;
                //currently does not correctly display if numAlleles > 2
                if (numAlleles == 0) {
                    sb.append("NA"); //if data does not exist
                }
                for (int i = 0; i < Math.min(2, numAlleles); i++) {
                    if (i > 0) {
                        sb.append('/');
                    }
                    sb.append(alleles[i]); //alleles
                }
                sb.append(delimChar);
                sb.append(alignment.getLocusName(site));
                sb.append(delimChar);
                sb.append(alignment.getPositionInLocus(site));
                sb.append(delimChar);
                sb.append("+"); //strand
                sb.append(delimChar);
                sb.append("NA"); //assembly# not supported
                sb.append(delimChar);
                sb.append("NA"); //center unavailable
                sb.append(delimChar);
                sb.append("NA"); //protLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //assayLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //panelLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //QCcode unavailable
                sb.append(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
                        String baseIUPAC = alignment.getBaseAsString(taxa, site);
                        sb.append(baseIUPAC);
                    } else {
                        String[] b = alignment.getBaseAsStringArray(taxa, site);
                        if (b.length == 1) {
                            sb.append(b[0]);
                            sb.append(b[0]);
                        } else {
                            sb.append(b[0]);
                            sb.append(b[1]);
                        }
                    }
                    if (taxa != (numTaxa - 1)) {
                        sb.append(delimChar);
                    }
                }
                bw.write(sb.toString());
                bw.write("\n");
            }
            bw.close();
        } catch (Exception e) {
            System.out.println("Error writing writeToHapmap: " + e);
            e.printStackTrace();
        }
    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error checking
     * @param alignment array of alignemnts
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static void writeToHapmap(Alignment alignment, AlignmentMask mask, boolean diploid, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".hmp.txt")), 1000000);
            bw.write("rs#");
            bw.write(delimChar);
            bw.write("alleles");
            bw.write(delimChar);
            bw.write("chrom");
            bw.write(delimChar);
            bw.write("pos");
            bw.write(delimChar);
            bw.write("strand");
            bw.write(delimChar);
            bw.write("assembly#");
            bw.write(delimChar);
            bw.write("center");
            bw.write(delimChar);
            bw.write("protLSID");
            bw.write(delimChar);
            bw.write("assayLSID");
            bw.write(delimChar);
            bw.write("panelLSID");
            bw.write(delimChar);
            bw.write("QCcode");
            bw.write(delimChar);
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {//finish filling out first row
                //not completely sure this does what I want, I need to access the
                //accession name from every alleleBLOB in bytes [52-201] but there
                //doesn't seem to be a method to access that in Alignment
                String sequenceID = alignment.getIdGroup().getIdentifier(taxa).getFullName().trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
                StringBuilder sb = new StringBuilder(alignment.getSNPID(site));
                sb.append(delimChar);
                byte[] alleles = alignment.getAlleles(site);
                int numAlleles = alleles.length;
                //currently does not correctly display if numAlleles > 2
                if (numAlleles == 0) {
                    sb.append("NA"); //if data does not exist
                }
                for (int i = 0; i < Math.min(2, numAlleles); i++) {
                    if (i > 0) {
                        sb.append('/');
                    }
                    sb.append(alleles[i]); //alleles
                }
                sb.append(delimChar);
                sb.append(alignment.getLocusName(site));
                sb.append(delimChar);
                sb.append(alignment.getPositionInLocus(site));
                sb.append(delimChar);
                sb.append("+"); //strand
                sb.append(delimChar);
                sb.append("NA"); //assembly# not supported
                sb.append(delimChar);
                sb.append("NA"); //center unavailable
                sb.append(delimChar);
                sb.append("NA"); //protLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //assayLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //panelLSID unavailable
                sb.append(delimChar);
                sb.append("NA"); //QCcode unavailable
                sb.append(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
                        String baseIUPAC = alignment.getBaseAsString(taxa, site);
                        if (mask.getMask(taxa, site) == 0x0) {
                            sb.append(baseIUPAC);
                        } else if (mask.getMask(taxa, site) == 0x1) {
                            sb.append(baseIUPAC.toLowerCase());
                        }
                    } else {
                        String[] b = alignment.getBaseAsStringArray(taxa, site);
                        if (b.length == 1) {
                            if (mask.getMask(taxa, site) == 0x0) {
                                sb.append(b[0]);
                                sb.append(b[0]);
                            } else if (mask.getMask(taxa, site) == 0x1) {
                                sb.append(b[0].toLowerCase());
                                sb.append(b[0].toLowerCase());
                            }
                        } else {
                            if (mask.getMask(taxa, site) == 0x0) {
                                sb.append(b[0]);
                                sb.append(b[1]);
                            } else if (mask.getMask(taxa, site) == 0x1) {
                                sb.append(b[0].toLowerCase());
                                sb.append(b[1].toLowerCase());
                            }
                        }
                    }
                    if (taxa != (numTaxa - 1)) {
                        sb.append(delimChar);
                    }
                }
                bw.write(sb.toString());
                bw.write("\n");
            }
            bw.close();
        } catch (Exception e) {
            System.out.println("Error writing writeToHapmap: " + e);
            e.printStackTrace();
        }
    }

    /**
     * Writes given set of alignments to a set of Plink files
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static void writeToPlink(Alignment alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        try {
            BufferedWriter MAPbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".plk.map")), 1000000);
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.getLocusName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(alignment.getSNPID(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write("NA"); // genetic distance unavailable
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.getPositionInLocus(site))); // position
                MAPbw.write("\n");
            }
            MAPbw.close();
            BufferedWriter PEDbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".plk.ped")), 1000000);
            // Compiled : Pattern
            Pattern splitter = Pattern.compile(":");
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] name = splitter.split(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim());
                if (name.length != 1) {
                    PEDbw.write(name[1]); // namelvl 1 if is available
                } else {
                    PEDbw.write("NA");
                }
                PEDbw.write(delimChar);
                PEDbw.write(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim()); // namelvl 0
                PEDbw.write(delimChar);
                PEDbw.write("NA"); // paternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("NA"); // maternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("other"); // gender is both
                PEDbw.write(delimChar);
                PEDbw.write("NA"); // phenotype unavailable, might have to change to "-9" for missing affection status
                PEDbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    String[] b = getSNPValueForPlink(alignment.getBaseAsStringArray(taxa, site));
                    PEDbw.write(b[0]);
                    PEDbw.write(delimChar);
                    PEDbw.write(b[b.length - 1]);
                    if (site != numSites - 1) {
                        PEDbw.write(delimChar);
                    }
                }
                PEDbw.write("\n");
            }
            PEDbw.close();
        } catch (Exception e) {
            myLogger.error("Error writing writeToPlink: " + e);
            throw new IllegalStateException("writeToPlink: " + e.getMessage());
        }
    }

    private static String[] getSNPValueForPlink(String[] base) {
        for (int i = 0; i < base.length; i++) {
            if (base[i].equals("N")) {
                base[i] = "0";
            } else if (base[i].equals("0")) {
                base[i] = "D";
            }
        }
        return base;
    }

    /**
     * Writes given set of alignments to a set of Flapjack files
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static void writeToFlapjack(Alignment alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        try {
            BufferedWriter MAPbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".flpjk.map")), 1000000);
            BufferedWriter DATbw = new BufferedWriter(new FileWriter(Utils.addSuffixIfNeeded(filename, ".flpjk.geno")), 1000000);
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.getSNPID(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write(alignment.getLocusName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.getPositionInLocus(site))); // position
                MAPbw.write("\n");
                DATbw.write(delimChar);
                DATbw.write(alignment.getSNPID(site));
            }
            MAPbw.close();
            DATbw.write("\n");
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                DATbw.write(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim());
                DATbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    String[] b = alignment.getBaseAsStringArray(taxa, site);
                    b = getSNPValueForFlapJack(b);
                    if (b.length == 1) {
                        DATbw.write(b[0]);
                    } else if (b.length == 2) {
                        DATbw.write(b[0]);
                        DATbw.write('/');
                        DATbw.write(b[1]);
                    } else if (b.length > 2) {
                        DATbw.write(b[0]);
                        DATbw.write('/');
                        DATbw.write(b[1]);
                        DATbw.write('/');
                        DATbw.write(b[2]);
                    }
                    if (site != numSites - 1) {
                        DATbw.write(delimChar);
                    }
                }
                DATbw.write("\n");
            }
            DATbw.close();
        } catch (Exception e) {
            System.out.println("Error writing writeToFlapjack: " + e);
        }
    }

    private static String[] getSNPValueForFlapJack(String[] base) {
        for (int i = 0; i < base.length; i++) {
            if (base[i].equals("N")) {
                base[i] = "-";
            }
        }
        return base;
    }

    public static void saveDelimitedAlignment(Alignment theAlignment, String delimit, String saveFile) {

        if ((saveFile == null) || (saveFile.length() == 0)) {
            return;
        }
        saveFile = Utils.addSuffixIfNeeded(saveFile, ".txt");
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(new File(saveFile));
            bw = new BufferedWriter(fw);

            bw.write("Taxa");
            int numSites = theAlignment.getSiteCount();
            for (int j = 0; j < numSites; j++) {
                bw.write(delimit);
                bw.write(String.valueOf(theAlignment.getPositionInLocus(j)));
            }
            bw.write("\n");

            for (int r = 0, n = theAlignment.getSequenceCount(); r < n; r++) {
                bw.write(theAlignment.getIdGroup().getIdentifier(r).getFullName());
                for (int i = 0; i < numSites; i++) {
                    bw.write(delimit);
                    bw.write(theAlignment.getBaseAsString(r, i));
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            System.out.println("AlignmentUtils: saveDelimitedAlignment: problem writing file: " + e.getMessage());
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    /** print alignment (in PHYLIP SEQUENTIAL format) */
    public static void printSequential(Alignment a, PrintWriter out) {
        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount() + "  S");

        // Print sequences
        for (int s = 0; s < a.getSequenceCount(); s++) {
            int n = 0;
            while (n < a.getSiteCount()) {
                if (n == 0) {
                    format.displayLabel(out,
                            a.getIdGroup().getIdentifier(s).getName(), 10);
                    out.print("     ");
                } else {
                    out.print("               ");
                }
                printNextSites(a, out, false, s, n, 50);
                out.println();
                n += 50;
            }
        }
    }

    /** print alignment (in PHYLIP 3.4 INTERLEAVED format) */
    public static void printInterleaved(Alignment a, PrintWriter out) {
        int n = 0;

        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount());

        // Print sequences
        while (n < a.getSiteCount()) {
            for (int s = 0; s < a.getSequenceCount(); s++) {
                if (n == 0) {
                    format.displayLabel(out,
                            a.getIdGroup().getIdentifier(s).getName(), 10);
                    out.print("     ");
                } else {
                    out.print("               ");
                }
                printNextSites(a, out, true, s, n, 50);
                out.println();
            }
            out.println();
            n += 50;
        }
    }

    /** Print alignment (in CLUSTAL W format) */
    public static void printCLUSTALW(Alignment a, PrintWriter out) {
        int n = 0;

        // CLUSTAL W header line
        out.println("CLUSTAL W multiple sequence alignment");
        out.println();

        // Print sequences
        while (n < a.getSiteCount()) {
            out.println();
            for (int s = 0; s < a.getSequenceCount(); s++) {
                format.displayLabel(out, a.getIdGroup().getIdentifier(s).getName(), 10);
                out.print("     ");

                printNextSites(a, out, false, s, n, 50);
                out.println();
            }
            // Blanks in status line are necessary for some parsers)
            out.println("               ");
            n += 50;
        }
    }

    private static void printNextSites(Alignment a, PrintWriter out, boolean chunked, int seq, int start, int num) {
        // Print next num characters
        for (int i = 0; (i < num) && (start + i < a.getSiteCount()); i++) {
            // Chunks of 10 characters
            if (i % 10 == 0 && i != 0 && chunked) {
                out.print(' ');
            }
            out.print(a.getBaseAsString(seq, start + i));
        }
    }
}
