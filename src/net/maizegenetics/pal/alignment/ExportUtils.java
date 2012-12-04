/*
 * ExportUtils
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;
import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 * The class exports PAL alignment data types to various file formats.
 *
 * @author Jon, Terry, Ed
 */
public class ExportUtils {

    private static final Logger myLogger = Logger.getLogger(ExportUtils.class);
    private static FormattedOutput format = FormattedOutput.getInstance();

    private ExportUtils() {
        // Utility Class - do not instantiate.
    }

    public static String writeToHDF5(Alignment a, String newHDF5file) {

        a = AlignmentUtils.optimizeForSites(a);
        a = AlignmentUtils.optimizeForTaxa(a);
        IHDF5Writer h5w = null;
        try {

            int numSites = a.getSiteCount();
            int numTaxa = a.getSequenceCount();

            newHDF5file = Utils.addSuffixIfNeeded(newHDF5file, ".hmp.h5");
            File hdf5File = new File(newHDF5file);
            if (hdf5File.exists()) {
                throw new IllegalArgumentException("ExportUtils: writeToHDF5: File already exists: " + newHDF5file);
            }
            IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
            myLogger.info("Writing HDF5 file: " + newHDF5file);
            config.overwrite();
            config.dontUseExtendableDataTypes();
            config.useUTF8CharacterEncoding();
            h5w = config.writer();
            
            h5w.setIntAttribute(HDF5Constants.DEFAULT_ATTRIBUTES_PATH, HDF5Constants.MAX_NUM_ALLELES, a.getMaxNumAlleles());
            
            h5w.setBooleanAttribute(HDF5Constants.DEFAULT_ATTRIBUTES_PATH, HDF5Constants.RETAIN_RARE_ALLELES, a.retainsRareAlleles());

            h5w.setIntAttribute(HDF5Constants.DEFAULT_ATTRIBUTES_PATH, HDF5Constants.NUM_TAXA, numTaxa);

            int numSBitWords = a.getAllelePresenceForAllTaxa(0, 0).getNumWords();
            h5w.setIntAttribute(HDF5Constants.DEFAULT_ATTRIBUTES_PATH, HDF5Constants.NUM_SBIT_WORDS, numSBitWords);

            int numTBitWords = a.getAllelePresenceForAllSites(0, 0).getNumWords();
            h5w.setIntAttribute(HDF5Constants.DEFAULT_ATTRIBUTES_PATH, HDF5Constants.NUM_TBIT_WORDS, numTBitWords);

            String[][] aEncodings = a.getAlleleEncodings();
            //int numEncodings = aEncodings.length;
            //numEncodings = 3;
            myLogger.info(Arrays.deepToString(aEncodings));
            h5w.writeStringArray(HDF5Constants.ALLELE_STATES, aEncodings[0]);
            //MDArray<String> alleleEncodings = new MDArray<String>(String.class, new int[]{numEncodings, 16});
            //String[] flatAlleleEncodings = new String[numEncodings * 16];
            //int count = 0;
            //for (int s = 0; s < numEncodings; s++) {
            //    for (int x = 0; x < 16; x++) {
                    //alleleEncodings.set(aEncodings[0][x], s, x);
                    //alleleEncodings.set(String.valueOf(s) + ":" + String.valueOf(x), s, x);
            //        System.out.println("flat: " + count + ": " + aEncodings[0][x]);
            //        flatAlleleEncodings[count++] = aEncodings[0][x];
            //    }
            //}
            //MDArray<String> alleleEncodings = new MDArray<String>(a.getAlleleEncodings(), new int[]{numSites, 16});
            //MDArray<String> alleleEncodings = new MDArray<String>(flatAlleleEncodings, new int[]{numEncodings, 16});
            //h5w.writeStringMDArray(HDF5Constants.ALLELE_STATES, alleleEncodings);
            
            h5w.writeStringArray(HDF5Constants.SNP_IDS, a.getSNPIDs());

            h5w.createGroup(HDF5Constants.SBIT);
            h5w.setIntAttribute(HDF5Constants.DEFAULT_ATTRIBUTES_PATH, HDF5Constants.NUM_SITES, numSites);

            String[] lociNames = new String[a.getNumLoci()];
            Locus[] loci = a.getLoci();
            for (int i = 0; i < a.getNumLoci(); i++) {
                lociNames[i] = loci[i].getName();
            }
            h5w.createStringVariableLengthArray(HDF5Constants.LOCI, a.getNumLoci());
            h5w.writeStringVariableLengthArray(HDF5Constants.LOCI, lociNames);

            h5w.createIntArray(HDF5Constants.LOCUS_OFFSETS, a.getNumLoci());
            h5w.writeIntArray(HDF5Constants.LOCUS_OFFSETS, a.getLociOffsets());

            h5w.createIntArray(HDF5Constants.POSITIONS, numSites);
            h5w.writeIntArray(HDF5Constants.POSITIONS, a.getPhysicalPositions());

            h5w.createByteMatrix(HDF5Constants.ALLELES, a.getSiteCount(), a.getMaxNumAlleles());
            byte[][] alleles = new byte[numSites][a.getMaxNumAlleles()];
            for (int i = 0; i < numSites; i++) {
                alleles[i] = a.getAlleles(i);
            }
            h5w.writeByteMatrix(HDF5Constants.ALLELES, alleles);

            String[] tn = new String[numTaxa];
            for (int i = 0; i < tn.length; i++) {
                tn[i] = a.getFullTaxaName(i);
            }
            h5w.createStringVariableLengthArray(HDF5Constants.TAXA, numTaxa);
            h5w.writeStringVariableLengthArray(HDF5Constants.TAXA, tn);

            for (int aNum = 0; aNum < a.getTotalNumAlleles(); aNum++) {

                String currentSBitPath = HDF5Constants.SBIT + "/" + aNum;
                h5w.createLongMatrix(currentSBitPath, numSites, numSBitWords, 1, numSBitWords);
                for (int i = 0; i < numSites; i++) {
                    long[][] lg = new long[1][numSBitWords];
                    lg[0] = a.getAllelePresenceForAllTaxa(i, aNum).getBits();
                    h5w.writeLongMatrixBlockWithOffset(currentSBitPath, lg, i, 0);
                }

                String currentTBitPath = HDF5Constants.TBIT + "/" + aNum;
                h5w.createLongMatrix(currentTBitPath, numTaxa, numTBitWords, 1, numTBitWords);
                for (int i = 0; i < numTaxa; i++) {
                    long[][] lg = new long[1][numTBitWords];
                    lg[0] = a.getAllelePresenceForAllSites(i, aNum).getBits();
                    h5w.writeLongMatrixBlockWithOffset(currentTBitPath, lg, i, 0);
                }

            }

            return newHDF5file;

        } finally {
            try {
                h5w.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error
     * checking
     *
     * @param alignment array of alignments
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static String writeToHapmap(Alignment alignment, boolean diploid, String filename, char delimChar, ProgressListener listener) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(filename, ".hmp.txt", new String[]{".hmp.txt", ".hmp.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
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
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                //finish filling out first row
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
                bw.write(alignment.getSNPID(site));
                bw.write(delimChar);
                byte[] alleles = alignment.getAlleles(site);
                int numAlleles = alleles.length;
                if (numAlleles == 0) {
                    bw.write("NA"); //if data does not exist
                } else if (numAlleles == 1) {
                    bw.write(alignment.getBaseAsString(site, alleles[0]));
                } else {
                    bw.write(alignment.getBaseAsString(site, alleles[0]));
                    if (alleles[1] != Alignment.UNKNOWN_ALLELE) {
                        bw.write('/');
                        bw.write(alignment.getBaseAsString(site, alleles[1]));
                    }
                }
                bw.write(delimChar);
                bw.write(alignment.getLocusName(site));
                bw.write(delimChar);
                bw.write(String.valueOf(alignment.getPositionInLocus(site)));
                bw.write(delimChar);
                bw.write("+"); //strand
                bw.write(delimChar);
                bw.write("NA"); //assembly# not supported
                bw.write(delimChar);
                bw.write("NA"); //center unavailable
                bw.write(delimChar);
                bw.write("NA"); //protLSID unavailable
                bw.write(delimChar);
                bw.write("NA"); //assayLSID unavailable
                bw.write(delimChar);
                bw.write("NA"); //panelLSID unavailable
                bw.write(delimChar);
                bw.write("NA"); //QCcode unavailable
                bw.write(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
                        String baseIUPAC = null;
                        try {
                            baseIUPAC = alignment.getBaseAsString(taxa, site);
                        } catch (Exception e) {
                            String[] b = alignment.getBaseAsStringArray(taxa, site);
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + "\nTry Exporting as Diploid Values.");
                        }
                        if (baseIUPAC == null) {
                            String[] b = alignment.getBaseAsStringArray(taxa, site);
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + "\nTry Exporting as Diploid Values.");
                        }
                        bw.write(baseIUPAC);
                    } else {
                        String[] b = alignment.getBaseAsStringArray(taxa, site);
                        if (b.length == 1) {
                            bw.write(b[0]);
                            bw.write(b[0]);
                        } else {
                            bw.write(b[0]);
                            bw.write(b[1]);
                        }
                    }
                    if (taxa != (numTaxa - 1)) {
                        bw.write(delimChar);
                    }
                }
                bw.write("\n");

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }
            return fullFileName;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing Hapmap file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error
     * checking
     *
     * @param alignment array of alignments
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static String writeToHapmap(Alignment alignment, AlignmentMask mask, boolean diploid, String filename, char delimChar, ProgressListener listener) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(filename, ".hmp.txt", new String[]{".hmp.txt", ".hmp.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
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
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                //finish filling out first row
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
                bw.write(alignment.getSNPID(site));
                bw.write(delimChar);
                byte[] alleles = alignment.getAlleles(site);
                int numAlleles = alleles.length;
                if (numAlleles == 0) {
                    bw.write("NA"); //if data does not exist
                } else if (numAlleles == 1) {
                    bw.write(alignment.getBaseAsString(site, alleles[0]));
                } else {
                    bw.write(alignment.getBaseAsString(site, alleles[0]));
                    if (alleles[1] != Alignment.UNKNOWN_ALLELE) {
                        bw.write('/');
                        bw.write(alignment.getBaseAsString(site, alleles[1]));
                    }
                }
                bw.write(delimChar);
                bw.write(alignment.getLocusName(site));
                bw.write(delimChar);
                bw.write(String.valueOf(alignment.getPositionInLocus(site)));
                bw.write(delimChar);
                bw.write("+"); //strand
                bw.write(delimChar);
                bw.write("NA"); //assembly# not supported
                bw.write(delimChar);
                bw.write("NA"); //center unavailable
                bw.write(delimChar);
                bw.write("NA"); //protLSID unavailable
                bw.write(delimChar);
                bw.write("NA"); //assayLSID unavailable
                bw.write(delimChar);
                bw.write("NA"); //panelLSID unavailable
                bw.write(delimChar);
                bw.write("NA"); //QCcode unavailable
                bw.write(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
                        String baseIUPAC = alignment.getBaseAsString(taxa, site);
                        if (mask.getMask(taxa, site) == 0x0) {
                            bw.write(baseIUPAC);
                        } else if (mask.getMask(taxa, site) == 0x1) {
                            bw.write(baseIUPAC.toLowerCase());
                        }
                    } else {
                        String[] b = alignment.getBaseAsStringArray(taxa, site);
                        if (b.length == 1) {
                            if (mask.getMask(taxa, site) == 0x0) {
                                bw.write(b[0]);
                                bw.write(b[0]);
                            } else if (mask.getMask(taxa, site) == 0x1) {
                                bw.write(b[0].toLowerCase());
                                bw.write(b[0].toLowerCase());
                            }
                        } else {
                            if (mask.getMask(taxa, site) == 0x0) {
                                bw.write(b[0]);
                                bw.write(b[1]);
                            } else if (mask.getMask(taxa, site) == 0x1) {
                                bw.write(b[0].toLowerCase());
                                bw.write(b[1].toLowerCase());
                            }
                        }
                    }
                    if (taxa != (numTaxa - 1)) {
                        bw.write(delimChar);
                    }
                }
                bw.write("\n");

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }
            return fullFileName;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing Hapmap file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Writes given set of alignments to a set of Plink files
     *
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static String writeToPlink(Alignment alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter MAPbw = null;
        BufferedWriter PEDbw = null;
        String mapFileName = Utils.addSuffixIfNeeded(filename, ".plk.map");
        String pedFileName = Utils.addSuffixIfNeeded(filename, ".plk.ped");
        try {
            MAPbw = new BufferedWriter(new FileWriter(mapFileName), 1000000);
            int numSites = alignment.getSiteCount();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.getLocusName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(alignment.getSNPID(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write("-9"); // genetic distance unavailable
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.getPositionInLocus(site))); // position
                MAPbw.write("\n");
            }
            MAPbw.close();

            PEDbw = new BufferedWriter(new FileWriter(pedFileName), 1000000);
            // Compiled : Pattern
            Pattern splitter = Pattern.compile(":");
            int numTaxa = alignment.getSequenceCount();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] name = splitter.split(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim());
                if (name.length != 1) {
                    PEDbw.write(name[1]); // namelvl 1 if is available
                } else {
                    PEDbw.write("-9");
                }
                PEDbw.write(delimChar);
                PEDbw.write(alignment.getIdGroup().getIdentifier(taxa).getFullName().trim()); // namelvl 0
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // paternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // maternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // gender is both
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // phenotype unavailable, might have to change to "-9" for missing affection status
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
            return mapFileName + " and " + pedFileName;
        } catch (Exception e) {
            myLogger.error("Error writing Plink files: " + mapFileName + " and " + pedFileName + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Plink files: " + mapFileName + " and " + pedFileName + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                PEDbw.close();
            } catch (Exception e) {
                // do nothing
            }
            try {
                MAPbw.close();
            } catch (Exception e) {
                // do nothing
            }
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
     *
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static String writeToFlapjack(Alignment alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        String mapFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.map");
        String genoFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.geno");
        try {

            BufferedWriter MAPbw = new BufferedWriter(new FileWriter(mapFileName), 1000000);
            BufferedWriter DATbw = new BufferedWriter(new FileWriter(genoFileName), 1000000);
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
            return mapFileName + " and " + genoFileName;

        } catch (Exception e) {
            myLogger.error("Error writing Flapjack files: " + mapFileName + " and " + genoFileName + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Flapjack files: " + mapFileName + " and " + genoFileName + ": " + ExceptionUtils.getExceptionCauses(e));
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

    public static String saveDelimitedAlignment(Alignment theAlignment, String delimit, String saveFile) {

        if ((saveFile == null) || (saveFile.length() == 0)) {
            return null;
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

            return saveFile;

        } catch (Exception e) {
            myLogger.error("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    /**
     * print alignment (in PHYLIP SEQUENTIAL format)
     */
    public static void printSequential(Alignment a, PrintWriter out) {
        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount() + "  S");

        // Print sequences
        for (int s = 0; s < a.getSequenceCount(); s++) {
            int n = 0;
            while (n < a.getSiteCount()) {
                if (n == 0) {
                    format.displayLabel(out, a.getIdGroup().getIdentifier(s).getName(), 10);
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

    /**
     * print alignment (in PHYLIP 3.4 INTERLEAVED format)
     */
    public static void printInterleaved(Alignment a, PrintWriter out) {
        int n = 0;

        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount());

        // Print sequences
        while (n < a.getSiteCount()) {
            for (int s = 0; s < a.getSequenceCount(); s++) {
                if (n == 0) {
                    format.displayLabel(out, a.getIdGroup().getIdentifier(s).getName(), 10);
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

    /**
     * Print alignment (in CLUSTAL W format)
     */
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

    public static String writeAlignmentToSerialGZ(Alignment sba, String outFile) {

        long time = System.currentTimeMillis();

        File theFile = null;
        FileOutputStream fos = null;
        GZIPOutputStream gz = null;
        ObjectOutputStream oos = null;
        try {
            theFile = new File(Utils.addSuffixIfNeeded(outFile, ".serial.gz"));
            fos = new FileOutputStream(theFile);
            gz = new GZIPOutputStream(fos);
            oos = new ObjectOutputStream(gz);
            oos.writeObject(sba);
            return theFile.getName();
        } catch (Exception e) {
            e.printStackTrace();
            myLogger.error("Error writing Serial GZ: " + theFile.getName() + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Serial GZ: " + theFile.getName() + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                oos.flush();
                oos.close();
                gz.close();
                fos.close();
            } catch (Exception e) {
                // do nothing
            }
            myLogger.info("writeAlignmentToSerialGZ: " + theFile.toString() + "  Time: " + (System.currentTimeMillis() - time));
        }

    }
}
