package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.File;

import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import org.apache.log4j.Logger;

/**
 * Basic filters needed for removing bad sites and taxa from GBS pipelines
 * @author edbuckler
 */
public class GBSHapMapFiltersPlugin extends AbstractPlugin {

    int startChromosome = 1, endChromosome = 10;
    private ArgsEngine myArgsEngine = null;
    private static final Logger myLogger = Logger.getLogger(GBSHapMapFiltersPlugin.class);
    private String suppliedInputFileName, suppliedOutputFileName, infile, outfile;
    private double minF = -2.0, minMAF = 0, maxMAF = 1, minPresence = 0;
    private boolean usePedigree = false;
    HashMap<String, Double> taxaFs = null;
    private boolean hLD = false;
    private double minR2 = 0.01;
    private double minBonP = 0.01;
    String[] lowCoverageTaxa = null;

    public GBSHapMapFiltersPlugin() {
        super(null, false);
    }

    public GBSHapMapFiltersPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        for (int chr = startChromosome; chr <= endChromosome; chr++) {
            infile = suppliedInputFileName.replace("+", "" + chr);
            outfile = suppliedOutputFileName.replace("+", "" + chr);
            myLogger.info("Reading: " + infile);
            Alignment a;
            try {
                a = ImportUtils.readFromHapmap(infile, null);
            } catch (Exception e) {
                myLogger.info("Could not read input hapmap file for chr" + chr + ":\n\t" + infile + "\n\tSkipping...");
                continue;
            }
            myLogger.info("Original Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount() == 0) {
                continue;
            }
            double realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, false);
            double randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, true, false);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a, false);
            myLogger.info("Ratio of RandomToReal:" + randomDist / realDist);
            if (myArgsEngine.getBoolean("-mnTCov")) {
                double tCov = Double.parseDouble(myArgsEngine.getString("-mnTCov"));
                if (lowCoverageTaxa == null) {
                    lowCoverageTaxa = getLowCoverageLines(a, tCov);
                }  // Note: lowCoverageTaxa is based upon the startChromosome only
                IdGroup keepTaxa = AlignmentFilterByGBSUtils.getFilteredIdGroupByName(a.getIdGroup(), lowCoverageTaxa, false);
                a = FilterAlignment.getInstance(a, keepTaxa);
                myLogger.info("TaxaFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
                if (a.getSiteCount() == 0) {
                    continue;
                }
            }
            realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, false);
            randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, true, false);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a, false);
            myLogger.info("Ratio of RandomToReal:" + randomDist / realDist);
            int minCount = (int) Math.round(a.getSequenceCount() * minPresence);
            if (usePedigree) {
                // filter the sites for minCount, minMAF and maxMAF (but not minF) based on all of the taxa
                int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, false, -2.0, minCount, minMAF, maxMAF);
                a = FilterAlignment.getInstance(a, goodLowHetSites);

                // filter the sites for minF only based only on the taxa with expectedF >= minF
                String[] highExpectedFTaxa = getHighExpectedFTaxa(a);
                IdGroup highExpectedFTaxaIDGroup = AlignmentFilterByGBSUtils.getFilteredIdGroupByName(a.getIdGroup(), highExpectedFTaxa, true);
                Alignment inbredGenos = FilterAlignment.getInstance(a, highExpectedFTaxaIDGroup);
                int[] goodLowFSites = AlignmentFilterByGBSUtils.getLowHetSNPs(inbredGenos, false, minF, 0, -0.1, 2.0);
                inbredGenos = null;
                System.gc();
                a = FilterAlignment.getInstance(a, goodLowFSites);
            } else {
                int[] goodLowHetSites = AlignmentFilterByGBSUtils.getLowHetSNPs(a, false, minF, minCount, minMAF, maxMAF);
                a = FilterAlignment.getInstance(a, goodLowHetSites);
            }
            myLogger.info("SiteFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
            if (a.getSiteCount() == 0) {
                continue;
            }
            realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, false, true);
            randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(a, true, true, false);
            System.out.printf("%d %d %g %g %g %n", a.getSiteCount(), minCount, minF, minMAF, randomDist / realDist);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(a, false);
            ExportUtils.writeToHapmap(a, false, outfile, '\t', null);
            myLogger.info("File written after basic filtering:" + outfile);

            if (hLD) {
                a = ImportUtils.readFromHapmap(outfile, null);
                int[] gs = AlignmentFilterByGBSUtils.getGoodSitesByLD(a, minR2, minBonP, 128, 100, 20, false);
                a = FilterAlignment.getInstance(a, gs);
                myLogger.info("LDFiltered Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
                if (a.getSiteCount() == 0) {
                    continue;
                }
                ExportUtils.writeToHapmap(a, false, outfile, '\t', null);
                myLogger.info("File written after basic & LD filtering:" + outfile);
            }
        }
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\n\n\nThe GBSHapMapFiltersPlugin accepts the following options:\n"
                + "-hmp     Input HapMap file; use a plus sign (+) as a wild card character to "
                + "           specify multiple chromosome numbers.\n"
                + "-o       Output HapMap file\n"
                + "-mnTCov  Minimum taxa coverage (default: no filter)\n"
                + "-mnSCov  Minimum presence (default: no filter)\n"
                + "-mnF     Minimum F (inbreeding coefficient) (default -2.0 = no filter)\n"
                + "-p       Pedigree file containing full sample names (or expected names after merging) & expected inbreeding\n"
                + "         coefficient (F) for each.  Only taxa with expected F >= mnF used to calculate F = 1-Ho/He.\n"
                + "         (default: use ALL taxa to calculate F)\n"
                + "-mnMAF   Minimum minor allele frequency (default: 0.0 = no filter)\n"
                + "-mxMAF   Maximum minor allele frequency (default: 1.0 = no filter)\n"
                + "-hLD     Filter for high LD\n"
                + "-mnR2    Minimum R-square value for the LD filter (default: " + minR2 + ")\n"
                + "-mnBonP  Minimum Bonferroni-corrected p-value for the LD filter (default: " + minBonP + ")\n"
                + "-sC      Start chromosome (default: 1).\n"
                + "-eC      End chromosome (default: 10).\n\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-hmp", "-hmpFile", true);
            myArgsEngine.add("-o", "--outFile", true);
            myArgsEngine.add("-mnTCov", "--minTaxaCov", true);
            myArgsEngine.add("-mnSCov", "--minSiteCov", true);
            myArgsEngine.add("-mnF", "--minFInbreeding", true);
            myArgsEngine.add("-p", "--pedigree-file", true);
            myArgsEngine.add("-mnMAF", "--minMinorAlleleFreq", true);
            myArgsEngine.add("-mxMAF", "--maxinorAlleleFreq", true);
            myArgsEngine.add("-hLD", "--highLD", false);
            myArgsEngine.add("-mnR2", "--minRSquare", true);
            myArgsEngine.add("-mnBonP", "--minBonferronPForLD", true);
            myArgsEngine.add("-sC", "--startChrom", true);
            myArgsEngine.add("-eC", "--endChrom", true);
        }

        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-sC")) {
            startChromosome = Integer.parseInt(myArgsEngine.getString("-sC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide a start chromosome.\n");
        }
        if (myArgsEngine.getBoolean("-eC")) {
            endChromosome = Integer.parseInt(myArgsEngine.getString("-eC"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please provide an end chromosome.\n");
        }
        if (myArgsEngine.getBoolean("-mnSCov")) {
            minPresence = Double.parseDouble(myArgsEngine.getString("-mnSCov"));
        }
        if (myArgsEngine.getBoolean("-mnF")) {
            minF = Double.parseDouble(myArgsEngine.getString("-mnF"));
        }
        if (myArgsEngine.getBoolean("-p")) {
            String pedigreeFileStr = myArgsEngine.getString("-p");
            File pedigreeFile = new File(pedigreeFileStr);
            if (!pedigreeFile.exists() || !pedigreeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the pedigree input file (-p option: " + pedigreeFileStr + ").");
            }
            taxaFs = TagsToSNPByAlignmentPlugin.readTaxaFsFromFile(pedigreeFile);
            if (taxaFs == null) {
                throw new IllegalArgumentException("Problem reading the pedigree file. Progam aborted.");
            }
            usePedigree = true;
        }
        if (myArgsEngine.getBoolean("-mnMAF")) {
            minMAF = Double.parseDouble(myArgsEngine.getString("-mnMAF"));
        }
        if (myArgsEngine.getBoolean("-mxMAF")) {
            maxMAF = Double.parseDouble(myArgsEngine.getString("-mxMAF"));
        }
        if (myArgsEngine.getBoolean("-hLD")) {
            hLD = true;
        }
        if (myArgsEngine.getBoolean("-mnR2")) {
            if (hLD) {
                minR2 = Double.parseDouble(myArgsEngine.getString("-mnR2"));
            } else {
                printUsage();
                throw new IllegalArgumentException("The -mnR2 option requires that the -hLD option is invoked\n");
            }
        }
        if (myArgsEngine.getBoolean("-mnBonP")) {
            if (hLD) {
                minBonP = Double.parseDouble(myArgsEngine.getString("-mnBonP"));
            } else {
                printUsage();
                throw new IllegalArgumentException("The -mnBonP option requires that the -hLD option is invoked\n");
            }
        }
        if (myArgsEngine.getBoolean("-hmp")) {
            suppliedInputFileName = myArgsEngine.getString("-hmp");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a HapMap file to filter.\n");
        }
        if (myArgsEngine.getBoolean("-o")) {
            suppliedOutputFileName = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file name.\n");
        }
    }

    public static void tagmapTextToBinary() {
        //String infile="c:/output9m.txt";
        //String outfile="c:/output_binary.dat";
        String infile = "/Users/edbuckler/SolexaAnal/GBS/test/14FCGBS.tg.txt";
        String outfile = "/Users/edbuckler/SolexaAnal/GBS/test/14FCGBS.tg.bin";
        TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(infile, false);
        theTOPM.sortTable(true);
        myLogger.info(theTOPM.getTagCount());
        theTOPM.writeBinaryFile(new File(outfile));
        theTOPM.writeTextFile(new File("/Users/edbuckler/SolexaAnal/GBS/test/14FCGBS.tg.sort.txt"));
    }

    public static void removeDuplicateTagsOnMap() {
//	String infile="c:/sample_tags.txt";
        String infile = "/Users/edbuckler/SolexaAnal/GBS/reftags/14FCGBS.tg.txt";
        String outfile = "/Users/edbuckler/SolexaAnal/GBS/test/14FCGBS.tg.ndup.txt";
        TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(infile, false);
        TagsOnPhysicalMap theTOPMNoDup = new TagsOnPhysicalMap(theTOPM, true);

        myLogger.info(theTOPMNoDup.getTagCount());
        theTOPMNoDup.writeTextFile(new File(outfile));
        theTOPMNoDup.writeBinaryFile(new File(outfile.replace(".txt", ".bin")));
    }

    public static String[] getLowCoverageLines(Alignment a, double pCoverage) {
        ArrayList<String> lowLines = new ArrayList<String>();
        for (int i = 0; i < a.getSequenceCount(); i++) {
            int covered = 0;
            for (int j = 0; j < a.getSiteCount(); j++) {
                if (a.getBase(i, j) != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    covered++;
                }
            }
            double propCovered = (double) covered / (double) a.getSiteCount();
//            myLogger.info(a.getTaxaName(i)+":"+propCovered);
            if (propCovered < pCoverage) {
                lowLines.add(a.getIdGroup().getIdentifier(i).getFullName());
            }
        }
        String[] lowL = lowLines.toArray(new String[0]);
        return lowL;
    }

    private String[] getHighExpectedFTaxa(Alignment a) {
        ArrayList<String> highFLines = new ArrayList<String>();
        int nInbredTaxa = 0;
        for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
            String fullTaxonName = a.getFullTaxaName(taxon);
            if (taxaFs.containsKey(fullTaxonName)) {
                if (taxaFs.get(fullTaxonName) >= minF) {
                    highFLines.add(fullTaxonName);
                    nInbredTaxa++;
                }
            }
        }
        myLogger.info(nInbredTaxa + " taxa with an Expected F >= the mnF of " + minF + " were found in the pedigree file (-p option)");
        String[] highF = highFLines.toArray(new String[0]);
        return highF;
    }

    public static void fuseHapMapWithMultipleChromosomes() {
        String[] infileHMP = {"/Users/edbuckler/SolexaAnal/SNP55K/SNP55K_maize282_AGPv1_origStrand_20100513.hmp.txt",
            "/Users/edbuckler/SolexaAnal/GBS/test/HiSeq282wLowHet_110215.hmp.txt"};
        byte[] infilePrefix = {(byte) 'A', (byte) 'G'};
        String outfileHMP = "/Users/edbuckler/SolexaAnal/GBS/test/HS55K_110215.hmp.txt";
        Alignment[] a = new Alignment[infileHMP.length];
        IdGroup commonIDs = null;
        int totalSites = 0;
        for (int i = 0; i < a.length; i++) {
            myLogger.info("Reading:" + infileHMP[i]);
            a[i] = ImportUtils.readFromHapmap(infileHMP[i], null);
            if (i == 0) {
                commonIDs = a[i].getIdGroup();
            } else {
                commonIDs = IdGroupUtils.getCommonIds(commonIDs, a[i].getIdGroup());
            }
            totalSites += a[i].getSiteCount();
        }

        //MutableNucleotideAlignment theMSA = new MutableNucleotideAlignment(commonIDs, totalSites, a[0].getLoci());
        MutableNucleotideAlignment theMSA = MutableNucleotideAlignment.getInstance(commonIDs, totalSites);
        int currSite = 0;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].getSiteCount(); j++) {
                theMSA.setLocusOfSite(currSite, a[i].getLocus(j));
                theMSA.setPositionOfSite(currSite, a[i].getPositionInLocus(j));
                //theMSA.setSitePrefix(currSite, infilePrefix[i]);
                //    theMSA.setStrandOfSite(currSite, (byte)(a[i].isPositiveStrand(j)==0));
                currSite++;
            }
        }
        for (int t = 0; t < theMSA.getSequenceCount(); t++) {
            currSite = 0;
            for (int i = 0; i < a.length; i++) {
                int aID = a[i].getIdGroup().whichIdNumber(theMSA.getTaxaName(t));
                for (int j = 0; j < a[i].getSiteCount(); j++) {
                    theMSA.setBase(t, currSite, a[i].getBase(aID, j));
                    currSite++;
                }
            }
        }
        theMSA.clean();
        //theMSA.sortSiteByPhysicalPosition();
        ExportUtils.writeToHapmap(theMSA, false, outfileHMP, '\t', null);
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
