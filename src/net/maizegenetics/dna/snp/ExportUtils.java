/*
 * ExportUtils
 */
package net.maizegenetics.dna.snp;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import net.maizegenetics.util.FormattedOutput;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import net.maizegenetics.util.VCFUtil;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

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

    public static String writeToHDF5(GenotypeTable a, String newHDF5file) {
        IHDF5Writer h5w = null;
        try {

            int numSites = a.numberOfSites();
            int numTaxa = a.numberOfTaxa();

            newHDF5file = Utils.addSuffixIfNeeded(newHDF5file, ".hmp.h5");
            File hdf5File = new File(newHDF5file);
            if (hdf5File.exists()) {
                throw new IllegalArgumentException("ExportUtils: writeToHDF5: File already exists: " + newHDF5file);
            }
            IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
            myLogger.info("Writing HDF5 file: " + newHDF5file);
            config.overwrite();
            config.dontUseExtendableDataTypes();
            h5w = config.writer();

            //h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.MAX_NUM_ALLELES, a.maxNumAlleles());

            h5w.setBooleanAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.RETAIN_RARE_ALLELES, a.retainsRareAlleles());

            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TAXA, numTaxa);

            int numSBitWords = a.allelePresenceForAllTaxa(0, 0).getNumWords();
            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SBIT_WORDS, numSBitWords);

            int numTBitWords = a.allelePresenceForAllSites(0, 0).getNumWords();
            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TBIT_WORDS, numTBitWords);

            String[][] aEncodings = a.alleleDefinitions();
            //myLogger.info(Arrays.deepToString(aEncodings));
            int numEncodings = aEncodings.length;
            int numStates = aEncodings[0].length;
            MDArray<String> alleleEncodings = new MDArray<String>(String.class, new int[]{numEncodings, numStates});
            for (int s = 0; s < numEncodings; s++) {
                for (int x = 0; x < numStates; x++) {
                    alleleEncodings.set(aEncodings[s][x], s, x);
                }
            }

            h5w.createStringMDArray(HapMapHDF5Constants.ALLELE_STATES, 100, new int[]{numEncodings, numStates});
            h5w.writeStringMDArray(HapMapHDF5Constants.ALLELE_STATES, alleleEncodings);
            MDArray<String> alleleEncodingReadAgain = h5w.readStringMDArray(HapMapHDF5Constants.ALLELE_STATES);
            if (alleleEncodings.equals(alleleEncodingReadAgain) == false) {
                throw new IllegalStateException("ExportUtils: writeToHDF5: Mismatch Allele States, expected '" + alleleEncodings + "', found '" + alleleEncodingReadAgain + "'!");
            }

            //h5w.writeStringArray(HapMapHDF5Constants.SNP_IDS, a.getSNPIDs());

            h5w.createGroup(HapMapHDF5Constants.SBIT);
            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SITES, numSites);

            String[] lociNames = new String[a.numChromosomes()];
            Chromosome[] loci = a.chromosomes();
            for (int i = 0; i < a.numChromosomes(); i++) {
                lociNames[i] = loci[i].getName();
            }
            h5w.createStringVariableLengthArray(HapMapHDF5Constants.LOCI, a.numChromosomes());
            h5w.writeStringVariableLengthArray(HapMapHDF5Constants.LOCI, lociNames);

            h5w.createIntArray(HapMapHDF5Constants.LOCUS_OFFSETS, a.numChromosomes());
            h5w.writeIntArray(HapMapHDF5Constants.LOCUS_OFFSETS, a.chromosomesOffsets());

            h5w.createIntArray(HapMapHDF5Constants.POSITIONS, numSites);
            h5w.writeIntArray(HapMapHDF5Constants.POSITIONS, a.physicalPositions());

            //h5w.createByteMatrix(HapMapHDF5Constants.ALLELES, a.numberOfSites(), a.maxNumAlleles());
            //byte[][] alleles = new byte[numSites][a.maxNumAlleles()];
            //for (int i = 0; i < numSites; i++) {
            //    alleles[i] = a.alleles(i);
            //}
            //h5w.writeByteMatrix(HapMapHDF5Constants.ALLELES, alleles);

            String[] tn = new String[numTaxa];
            for (int i = 0; i < tn.length; i++) {
                tn[i] = a.taxaName(i);
            }
            h5w.createStringVariableLengthArray(HapMapHDF5Constants.TAXA, numTaxa);
            h5w.writeStringVariableLengthArray(HapMapHDF5Constants.TAXA, tn);

            //for (int aNum = 0; aNum < a.getTotalNumAlleles(); aNum++) {
            for (int aNum = 0; aNum < 6; aNum++) {

                String currentSBitPath = HapMapHDF5Constants.SBIT + "/" + aNum;
                h5w.createLongMatrix(currentSBitPath, numSites, numSBitWords, 1, numSBitWords);
                long[][] lgarray = new long[numSites][numSBitWords];
                for (int i = 0; i < numSites; i++) {
                    lgarray[i] = a.allelePresenceForAllTaxa(i, aNum).getBits();
                }
                h5w.writeLongMatrix(currentSBitPath, lgarray);

                String currentTBitPath = HapMapHDF5Constants.TBIT + "/" + aNum;
                h5w.createLongMatrix(currentTBitPath, numTaxa, numTBitWords, 1, numTBitWords);
                lgarray = new long[numTaxa][numTBitWords];
                for (int i = 0; i < numTaxa; i++) {
                    lgarray[i] = a.allelePresenceForAllSites(i, aNum).getBits();
                }
                h5w.writeLongMatrix(currentTBitPath, lgarray);

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
    
    public static String writeToMutableHDF5(GenotypeTable a, String newHDF5file) {
        return writeToMutableHDF5(a, newHDF5file, null, false);
    }
   
   /**
    * Exports a alignment into the Byte HDF5 format.  
    * @param a alignment to be exported
    * @param newHDF5file filename for the new file (should end with "hmp.h5")
    * @param exportTaxa  subset of taxa (if null exports ALL taxa)
    * @return 
    */ 
   public static String writeToMutableHDF5(GenotypeTable a, String newHDF5file, TaxaList exportTaxa, boolean keepDepth) {
        AlignmentBuilder aB=AlignmentBuilder.getTaxaIncremental(a.positions(),newHDF5file);
        if((exportTaxa!=null)&&(exportTaxa.numberOfTaxa()==0)) {aB.build(); return newHDF5file;}
        for (int t = 0; t < a.numberOfTaxa(); t++) {
              if((exportTaxa!=null)&&(!exportTaxa.contains(a.taxa().get(t)))) continue;  //taxon not in export list
              byte[] bases = a.genotypeAllSites(t);
              if (keepDepth==false) aB.addTaxon(new Taxon(a.taxaName(t)), bases, null);
              else {
                  //todo restore depth save
//                  MutableNucleotideAlignmentHDF5 m= (MutableNucleotideAlignmentHDF5) a;
//                  addA.addTaxon(new Taxon(a.taxaName(t)), bases, m.depthForAlleles(t));
              }
        }
        aB.build();
        return newHDF5file;
    }


    //TODO TAS-55 on new filters for HDF5 covers this.
//   public static String writeToMutableHDF5(Alignment a, String newHDF5file, int[] snpIndex) {
//        return writeToMutableHDF5(a, newHDF5file, snpIndex, null, false);
//   }
//
//   public static String writeToMutableHDF5(Alignment a, String newHDF5file, int[] snpIndex, TaxaList exportTaxa, boolean keepDepth) {
//        IHDF5Writer h5w = null;
//        try {
//            int numSites = snpIndex.length;
//            int numTaxa = a.numberOfTaxa();
//            newHDF5file = Utils.addSuffixIfNeeded(newHDF5file, "hmp.h5");
//            File hdf5File = new File(newHDF5file);
//            if (hdf5File.exists()) {
//                throw new IllegalArgumentException("ExportUtils: writeToMutableHDF5: File already exists: " + newHDF5file);
//            }
//            IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
//            myLogger.info("Writing Mutable HDF5 file: " + newHDF5file);
//            config.overwrite();
//            config.dontUseExtendableDataTypes();
//            h5w = config.writer();
//            HDF5IntStorageFeatures features = HDF5IntStorageFeatures.createDeflation(HDF5IntStorageFeatures.MAX_DEFLATION_LEVEL);
//            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.MAX_NUM_ALLELES, a.maxNumAlleles());
//            h5w.setBooleanAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.RETAIN_RARE_ALLELES, a.retainsRareAlleles());
//            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TAXA, numTaxa);
//            String[][] aEncodings = a.alleleDefinitions();
//            int numEncodings = aEncodings.length;
//            int numStates = aEncodings[0].length;
//            MDArray<String> alleleEncodings = new MDArray<String>(String.class, new int[]{numEncodings, numStates});
//            for (int s = 0; s < numEncodings; s++) {
//                for (int x = 0; x < numStates; x++) {
//                    alleleEncodings.set(aEncodings[s][x], s, x);
//                }
//            }
//            h5w.createStringMDArray(HapMapHDF5Constants.ALLELE_STATES, 100, new int[]{numEncodings, numStates});
//            h5w.writeStringMDArray(HapMapHDF5Constants.ALLELE_STATES, alleleEncodings);
//            MDArray<String> alleleEncodingReadAgain = h5w.readStringMDArray(HapMapHDF5Constants.ALLELE_STATES);
//            if (alleleEncodings.equals(alleleEncodingReadAgain) == false) {
//                throw new IllegalStateException("ExportUtils: writeToMutableHDF5: Mismatch Allele States, expected '" + alleleEncodings + "', found '" + alleleEncodingReadAgain + "'!");
//            }
//            String[] originalSNPID = a.getSNPIDs();
//            String[] outputSNPID = new String[snpIndex.length];
//            for (int i = 0; i < outputSNPID.length; i++) {
//                outputSNPID[i] = originalSNPID[snpIndex[i]];
//            }
//            h5w.writeStringArray(HapMapHDF5Constants.SNP_IDS, outputSNPID);  //consider adding compression
//
//            h5w.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SITES, numSites);
//
//            TreeSet<Chromosome> lset = new TreeSet();
//            for (int i = 0; i < numSites; i++) {
//                lset.add(a.chromosome(snpIndex[i]));
//            }
//            Chromosome[] outLoci = lset.toArray(new Chromosome[lset.size()]);
//            String[] lociNames = new String[outLoci.length];
//            HashMap<Chromosome, Integer> locusToIndex=new HashMap<Chromosome, Integer>(10);
//            for (int i = 0; i < lociNames.length; i++) {
//                lociNames[i] = outLoci[i].getName();
//                locusToIndex.put(outLoci[i],i);
//            }
//
//            h5w.createStringVariableLengthArray(HapMapHDF5Constants.LOCI, lociNames.length);
//            h5w.writeStringVariableLengthArray(HapMapHDF5Constants.LOCI, lociNames);
//
////            h5w.createIntArray(HapMapHDF5Constants.LOCUS_OFFSETS, a.getNumLoci());
////            h5w.writeIntArray(HapMapHDF5Constants.LOCUS_OFFSETS, a.getLociOffsets());
//
//            int[] locusIndicesArray = new int[snpIndex.length];
//            for (int i = 0; i < locusIndicesArray.length; i++) {
//                locusIndicesArray[i] = locusToIndex.get(a.chromosome(snpIndex[i]));
//            }
//
//            h5w.createIntArray(HapMapHDF5Constants.LOCUS_INDICES, a.numberOfSites(),features);
//            h5w.writeIntArray(HapMapHDF5Constants.LOCUS_INDICES, locusIndicesArray,features);
//
//            int[] originalPositions = a.physicalPositions();
//            int[] outPositions = new int[snpIndex.length];
//            for (int i = 0; i < snpIndex.length; i++) {
//                outPositions[i] = originalPositions[snpIndex[i]];
//            }
//            h5w.createIntArray(HapMapHDF5Constants.POSITIONS, numSites);
//            h5w.writeIntArray(HapMapHDF5Constants.POSITIONS, outPositions);
//
//            // Write Bases
//
//  //        HDF5IntStorageFeatures features = HDF5IntStorageFeatures.createDeflation(HDF5IntStorageFeatures.NO_DEFLATION_LEVEL);
//
//            h5w.createGroup(HapMapHDF5Constants.GENOTYPES);
//            if((exportTaxa!=null)&&(exportTaxa.numberOfTaxa()==0)) {h5w.close(); return newHDF5file;}
//            h5w.close();
//            MutableNucleotideAlignmentHDF5 addA=MutableNucleotideAlignmentHDF5.getInstance(newHDF5file);
//            for (int t = 0; t < numTaxa; t++) {
//                  if((exportTaxa!=null)&&(exportTaxa.whichIdNumber(a.taxaName(t))<0)) continue;  //taxon not in export list
//                  byte[] originalBases = a.genotypeAllSites(t);
//                  byte[] bases = new byte[snpIndex.length];
//                  for (int i = 0; i < bases.length; i++) {
//                      bases[i] = originalBases[snpIndex[i]];
//                  }
//
//                  if (keepDepth==false) addA.addTaxon(new Taxon(a.taxaName(t)), bases, null);
//                  else {
//                      MutableNucleotideAlignmentHDF5 m= (MutableNucleotideAlignmentHDF5) a;
//                      byte[][] originalDepth = m.depthForAlleles(t);
//                      byte[][] depth = new byte[originalDepth.length][snpIndex.length];
//                      for (int i = 0; i < depth.length; i++) {
//                          for (int j = 0; j < snpIndex.length; j++) {
//                              depth[i][j] = originalDepth[i][snpIndex[j]];
//                          }
//                      }
//                      addA.addTaxon(new Taxon(a.taxaName(t)), bases, depth);
//                  }
//            }
//            addA.clean();
//            return newHDF5file;
//
//        } finally {
//            try {
//                h5w.close();
//            } catch (Exception e) {
//                // do nothing
//            }
//        }
//    }

    
//         /**
//     * This merge multiple alignment together into one ByteNucleotideHDF5 File.
//     * This is designed for putting multiple chromosomes together into one whole genome file.  The
//     * sites must be exactly the same.
//     * <p>
//     * TODO:  Need to add option for merging depth.
//     * @param sourceFiles array of sourcefile names in ByteNucleotideHDF5 format with the new taxa
//     * @param targetFile name of ByteNucleotideHDF5 file to receive the new taxa
//     */
//    public static void addTaxaFromExistingByteHDF5File(String[] sourceFiles, String targetFile, boolean addDepth) {
//        System.out.println("Opening:"+targetFile);
//        MutableNucleotideAlignmentHDF5 trgA=(MutableNucleotideAlignmentHDF5)ImportUtils.readGuessFormat(targetFile, false);
//        for (String srcFile : sourceFiles) {
//           System.out.println("Opening:"+srcFile);
//            MutableNucleotideAlignmentHDF5 srcA=(MutableNucleotideAlignmentHDF5)ImportUtils.readGuessFormat(srcFile, false);
//            System.out.println("Target now has taxon:"+trgA.numberOfTaxa());
//            if(srcA.numberOfSites()!=trgA.numberOfSites()) {
//                throw new IllegalStateException("ExportUtils: addTaxaFromExistingByteHDF5File: Mismatch in number of sites");
//            }
//            System.out.println("Copying first taxon:"+srcA.taxaName(0));
//            for (int i = 0; i < srcA.numberOfTaxa(); i++) {
//                if (addDepth==true) trgA.addTaxon(srcA.taxa().getTaxon(i), srcA.genotypeAllSites(i), srcA.depthForAlleles(i));
//                else trgA.addTaxon(srcA.taxa().getTaxon(i), srcA.genotypeAllSites(i), null);
//            }
//        }
//
//        trgA.clean();
//        System.out.println("Target now has taxon:"+trgA.numberOfTaxa());
//    }

    /**
     * Writes multiple alignments to single Hapmap file. Currently no error
     * checking
     *
     * @param alignment array of alignments
     * @param diploid
     * @param filename
     * @param delimChar
     */
    public static String writeToHapmap(GenotypeTable alignment, boolean diploid, String filename, char delimChar, ProgressListener listener) {
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
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                //finish filling out first row
                //not completely sure this does what I want, I need to access the
                //accession name from every alleleBLOB in bytes [52-201] but there
                //doesn't seem to be a method to access that in Alignment
                String sequenceID = alignment.taxaName(taxa).trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                bw.write(alignment.siteName(site));
                bw.write(delimChar);
//                byte[] alleles = alignment.alleles(site); // doesn't work right for MutableVCFAlignment (always returns 3 alleles, even if no data)
//                int numAlleles = alleles.length;
                int[][] sortedAlleles = alignment.allelesSortedByFrequency(site); // which alleles are actually present among the genotypes
                int numAlleles = sortedAlleles[0].length;
                if (numAlleles == 0) {
                    bw.write("NA"); //if data does not exist
                } else if (numAlleles == 1) {
                    bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                } else {
                    bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                    for (int allele = 1; allele < sortedAlleles[0].length; allele++) {
                        if (sortedAlleles[0][allele] != GenotypeTable.UNKNOWN_ALLELE) {
                            bw.write('/');
                            bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][allele]));  // will write out a third allele if it exists
                        }
                    }
                }
                bw.write(delimChar);
                bw.write(alignment.chromosomeName(site));
                bw.write(delimChar);
                bw.write(String.valueOf(alignment.chromosomalPosition(site)));
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
                            baseIUPAC = alignment.genotypeAsString(taxa, site);
                        } catch (Exception e) {
                            String[] b = alignment.genotypeAsStringArray(taxa, site);
                            bw.close();
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                        }
                        if ((baseIUPAC == null) || baseIUPAC.equals("?")) {
                            String[] b = alignment.genotypeAsStringArray(taxa, site);
                            bw.close();
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                        }
                        bw.write(baseIUPAC);
                    } else {
                        String[] b = alignment.genotypeAsStringArray(taxa, site);
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
     * Writes given alignment to a VCF file
     *
     * @param alignment
     * @param filename
     * @return
     */
    public static String writeToVCF(GenotypeTable alignment, String filename, char delimChar) {
        //todo restore depth
        boolean hasDepth=false;  //in future test for this
        try {

            filename = Utils.addSuffixIfNeeded(filename, ".vcf", new String[]{".vcf", ".vcf.gz"});
            BufferedWriter bw = Utils.getBufferedWriter(filename);
            bw.write("##fileformat=VCFv4.0");
            bw.newLine();
            if (alignment.referenceGenotype(0) == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                bw.write("##Reference allele is not known. The major allele was used as reference allele.");
                bw.newLine();
            }
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
            bw.newLine();
            bw.write("####FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">");
            bw.newLine();
            bw.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
            bw.newLine();
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
            bw.newLine();


            bw.write("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">");
            bw.newLine();
            bw.write("#CHROM" + delimChar + "POS" + delimChar + "ID" + delimChar + "REF" + delimChar + "ALT" + delimChar + "QUAL" + delimChar + "FILTER" + delimChar + "INFO" + delimChar + "FORMAT");
            boolean refTaxon = false;
            for (int taxa = 0; taxa < alignment.numberOfTaxa(); taxa++) {
                String taxonName = alignment.taxaName(taxa).trim();
                if (taxa == 0 && taxonName.contentEquals("REFERENCE_GENOME")) {
                    refTaxon = true;
                } else {
                    bw.write(delimChar + taxonName);;
                }
            }
            bw.newLine();

            for (int site = 0; site < alignment.numberOfSites(); site++) {
                int[][] sortedAlleles = alignment.allelesSortedByFrequency(site); // which alleles are actually present among the genotypes


                int nAlleles = sortedAlleles[0].length;


                if (nAlleles == 0) {                                                  //used to be ==0
                    System.out.println("no alleles at: " + site + " " + alignment.chromosomalPosition(site));
                    continue;
                }

                byte refGeno = alignment.referenceGenotype(site);
                if (refGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    String myMajorAllele = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0][sortedAlleles[0][0]];
                    String MajorGenotype = myMajorAllele + myMajorAllele;
                    refGeno = NucleotideAlignmentConstants.getNucleotideDiploidByte(MajorGenotype);
                }
                byte refAllele = (byte) (refGeno & 0xF);  // converts from diploid to haploid allele (2nd allele)
                //System.out.println(alignment.chromosomalPosition(site) + " " + refAllele);
                byte[] alleleValues = null;
                if (hasDepth) {
                    alleleValues = alignment.allelesBySortType(GenotypeTable.ALLELE_SORT_TYPE.Depth, site); // storage order of the alleles in the alignment (myCommonAlleles & myAlleleDepth) (length always 3, EVEN IF THERE ARE ONLY 2 in the genos)
                } else {
                    alleleValues = alignment.allelesBySortType(GenotypeTable.ALLELE_SORT_TYPE.Frequency, site);
                    //if (nAlleles > alignment.maxNumAlleles()) {
                    //    nAlleles = alignment.maxNumAlleles();
                    //}
                }
                int[] alleleRedirect = new int[nAlleles]; // holds the indices of alleleValues in ref, alt1, [alt2] order (max 3 alleles)
                String refAlleleStr;
                int refUnknownOffset = 0;
                if (refGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {  // reference allele unknown - report the alleles in maj, min1, [min2] order
                    refUnknownOffset = 1;
                    refAlleleStr = ".";
                    int aRedirectIndex = 0;
                    for (int sortedAIndex = 0; sortedAIndex < nAlleles; sortedAIndex++) {
                        for (int aValuesIndex = 0; aValuesIndex < alleleValues.length; aValuesIndex++) {
                            if (alleleValues[aValuesIndex] == sortedAlleles[0][sortedAIndex]) {
                                alleleRedirect[aRedirectIndex] = aValuesIndex;
                                ++aRedirectIndex;
                                break;
                            }
                        }
                    }
                } else {  // refAllele known
                    refAlleleStr = NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.get(refGeno);
                    //refAlleleStr = String.valueOf(refGeno);

                    // check if the reference allele is found among the genotypes
                    boolean refAlleleAmongGenos = false;
                    boolean refAlleleInAllelValues = false;
                    for (int sortedAIndex = 0; sortedAIndex < nAlleles; sortedAIndex++) {
                        if (sortedAlleles[0][sortedAIndex] == refAllele) {
                            refAlleleAmongGenos = true;
                            for (int aValuesIndex = 0; aValuesIndex < alleleValues.length; aValuesIndex++) {
                                if (alleleValues[aValuesIndex] == refAllele) {
                                    refAlleleInAllelValues = true;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    if (refAlleleInAllelValues) {
                        // the refAllele index in alleleValues should be stored in alleleRedirect[0], and the remaining in desc order of freq
                        int aRedirectIndex = 1;
                        for (int sortedAIndex = 0; sortedAIndex < nAlleles; sortedAIndex++) {
                            for (int aValuesIndex = 0; aValuesIndex < alleleValues.length; aValuesIndex++) {
                                if (alleleValues[aValuesIndex] == sortedAlleles[0][sortedAIndex]) {
                                    if (alleleValues[aValuesIndex] == refAllele) {
                                        alleleRedirect[0] = aValuesIndex;
                                    } else {
                                        alleleRedirect[aRedirectIndex] = aValuesIndex;
                                        ++aRedirectIndex;
                                    }
                                    break;
                                }
                            }
                        }
                    } else {
                        // alleleRedirect[0] set to -1, the remaining in desc order of freq: maj, min1, [min2]
                        if (!refAlleleAmongGenos) {
                            alleleRedirect = new int[alleleRedirect.length + 1];
                        }
                        alleleRedirect[0] = -1;
                        int aRedirectIndex = 1;
                        for (int sortedAIndex = 0; sortedAIndex < nAlleles; sortedAIndex++) {
                            for (int aValuesIndex = 0; aValuesIndex < alleleValues.length; aValuesIndex++) {
                                if (alleleValues[aValuesIndex] == sortedAlleles[0][sortedAIndex]) {
                                    alleleRedirect[aRedirectIndex] = aValuesIndex;
                                    ++aRedirectIndex;
                                    break;
                                }
                            }
                        }
                    }
                }
                bw.write(alignment.chromosomeName(site)); // chromosome
                bw.write(delimChar);
                bw.write(alignment.chromosomalPosition(site) + ""); // position
                bw.write(delimChar);
                bw.write(alignment.siteName(site)); // site name
                bw.write(delimChar);
                bw.write(refAlleleStr); // ref allele
                bw.write(delimChar);

                StringBuilder altAllelesBuilder = new StringBuilder("");
                int firstAltAllele = (refAlleleStr == ".") ? 0 : 1;  // if ref allele is unknown, ALL alleles are alt
                for (int i = firstAltAllele; i < alleleRedirect.length; i++) {
                    if (i != firstAltAllele) {
                        altAllelesBuilder.append(",");
                    }
                    byte diploidByte = (byte) ((alleleValues[alleleRedirect[i]] << 4) | alleleValues[alleleRedirect[i]]);
                    String AlleleStr = NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.get(diploidByte);
                    if (AlleleStr == null) {
                        altAllelesBuilder.append(".");
                    } else {
                        altAllelesBuilder.append(AlleleStr);
                    }
                }
                String altAlleles = altAllelesBuilder.toString();
                altAlleles = (altAlleles.equals("")) ? "." : altAlleles;
                bw.write(altAlleles); // alt alleles
                bw.write(delimChar);

                bw.write("."); // qual score
                bw.write(delimChar);

                bw.write("PASS"); // filter
                bw.write(delimChar);

                if (hasDepth) {
                    int totalDepth = 0;
                    for (int i = 0; i < alignment.numberOfTaxa(); i++) {
                        byte[] depth = alignment.depthForAlleles(i, site);
                        for (int k = 0; k < depth.length; k++) {
                            if (depth[k] != -1) {
                                totalDepth += depth[k];
                            }
                        }
                    }
                    bw.write("DP=" + totalDepth); // DP
                } else {
                    bw.write("."); // DP
                }
                bw.write(delimChar);

                if (hasDepth) {
                    bw.write("GT:AD:DP:GQ:PL");
                } else {
                    bw.write("GT");
                }
                for (int taxa = 0; taxa < alignment.numberOfTaxa(); taxa++) {
                    if (taxa == 0 && refTaxon) {
                        continue;  // don't include REFERENCE_GENOME in vcf output
                    }
                    bw.write(delimChar);

                    // GT = genotype
                    String GTstr = "";
                    byte[] values = alignment.genotypeArray(taxa, site);

                    boolean genoOne = false;
                    if (values[0] == GenotypeTable.UNKNOWN_ALLELE) {
                        GTstr += "./";
                        genoOne = true;
                    } else {
                        for (int i = 0; i < alleleRedirect.length; i++) { // alleleRedirect stores the alleles in ref/alt1/[alt2] order (if no alt2,length=2)
                            if (i == 0 && alleleRedirect[i] == -1) {  // refAllele known but either not among genos or not in alleleValues
                                if (values[0] == refAllele) {

                                    GTstr += (i + refUnknownOffset) + "/";
                                    genoOne = true;
                                    break;
                                }
                            } else if (values[0] == alleleValues[alleleRedirect[i]]) {
                                GTstr += (i + refUnknownOffset) + "/";
                                genoOne = true;
                                break;
                            }
                        }
                    }
//                    if (!genoOne) { 
                    //bw.write("./.");
//                        if (values[0] == NucleotideAlignmentConstants.A_ALLELE) {
//                            bw.write("A/");
//                        } else if (values[0] == NucleotideAlignmentConstants.C_ALLELE) {
//                            bw.write("C/");
//                        } else if (values[0] == NucleotideAlignmentConstants.G_ALLELE) {
//                            bw.write("G/");
//                        } else if (values[0] == NucleotideAlignmentConstants.T_ALLELE) {
//                            bw.write("T/");
//                        } else if (values[0] == NucleotideAlignmentConstants.GAP_ALLELE) {
//                            bw.write("-/");
//                        } else {
//                            bw.write(values[0] + "/");
//                            // throw new IllegalArgumentException("Unknown allele value: " + alleleValues[i]);
//                        }
//                    }

                    boolean genoTwo = false;
                    if (values[1] == GenotypeTable.UNKNOWN_ALLELE) {
                        GTstr += ".";
                        genoTwo = true;
                    } else {
                        for (int i = 0; i < alleleRedirect.length; i++) { // alleleRedirect stores the alleles in ref/alt1/alt2 order (if no alt2,length=2)
                            if (i == 0 && alleleRedirect[i] == -1) {  // refAllele known but either not among genos or not in alleleValues
                                if (values[1] == refAllele) {
                                    GTstr += (i + refUnknownOffset) + "";
                                    genoTwo = true;
                                    break;
                                }
                            } else if (values[1] == alleleValues[alleleRedirect[i]]) {
                                GTstr += (i + refUnknownOffset) + "";
                                genoTwo = true;
                                break;
                            }
                        }
                    }
//                    if (!genoTwo) {
//                        if (values[1] == NucleotideAlignmentConstants.A_ALLELE) {
//                            bw.write("A");
//                        } else if (values[1] == NucleotideAlignmentConstants.C_ALLELE) {
//                            bw.write("C");
//                        } else if (values[1] == NucleotideAlignmentConstants.G_ALLELE) {
//                            bw.write("G");
//                        } else if (values[1] == NucleotideAlignmentConstants.T_ALLELE) {
//                            bw.write("T");
//                        } else if (values[1] == NucleotideAlignmentConstants.GAP_ALLELE) {
//                            bw.write("-");
//                        } else {
//                            bw.write(values[1] + "");
//                            // throw new IllegalArgumentException("Unknown allele value: " + alleleValues[i]);
//                        }
//                    }
                    if (genoOne && genoTwo) {
                        bw.write(GTstr);
                    } else {
                        bw.write("./.");
                    }
                    if (!(hasDepth)) {
                        continue;
                    }
                    bw.write(":");

                    // AD
                    byte[] siteAlleleDepths = alignment.depthForAlleles(taxa, site);

                    int siteTotalDepth = 0;
                    if (siteAlleleDepths.length != 0) {
                        for (int a = 0; a < alleleRedirect.length; a++) {
                            if (a == 0 && alleleRedirect[a] == -1) {  // refAllele known but is not among the genotypes
                                bw.write("0");
                                continue;
                            }
                            if (a != 0) {
                                bw.write(",");
                            }
                            bw.write(((int) (siteAlleleDepths[alleleRedirect[a]])) + "");
                            siteTotalDepth += siteAlleleDepths[alleleRedirect[a]];
                        }
                    } else {
                        bw.write(".,.,.");
                    }
                    bw.write(":");

                    // DP
                    bw.write(siteTotalDepth + "");
                    bw.write(":");

                    if (siteAlleleDepths.length != 0) {
                        int[] scores;
                        if (siteAlleleDepths.length == 1) {
                            int dep1 = siteAlleleDepths[alleleRedirect[0]] > 127 ? 127 : siteAlleleDepths[alleleRedirect[0]];
                            scores = VCFUtil.getScore(dep1, 0);
                        } else {
                            if (alleleRedirect[0] == -1) {
                                int dep1 = 0;
                                int dep2 = 0;
                                if (alleleRedirect.length > 2) {
                                    dep1 = siteAlleleDepths[alleleRedirect[1]] > 127 ? 127 : siteAlleleDepths[alleleRedirect[1]];
                                    dep2 = siteAlleleDepths[alleleRedirect[2]] > 127 ? 127 : siteAlleleDepths[alleleRedirect[2]];
                                } else if (alleleRedirect.length == 2) {
                                    dep1 = 0;
                                    dep2 = siteAlleleDepths[alleleRedirect[1]] > 127 ? 127 : siteAlleleDepths[alleleRedirect[1]];
                                }
                                scores = VCFUtil.getScore(dep1, dep2);
                            } else {
                                int dep1 = siteAlleleDepths[alleleRedirect[0]] > 127 ? 127 : siteAlleleDepths[alleleRedirect[0]];
                                int dep2 = 0;
                                if (alleleRedirect.length > 1) {
                                    dep2 = siteAlleleDepths[alleleRedirect[1]] > 127 ? 127 : siteAlleleDepths[alleleRedirect[1]];
                                }
                                scores = VCFUtil.getScore(dep1, dep2);
                            }
                        }

                        // GQ
                        if (scores == null) {
                            scores = new int[]{-1, -1, -1, -1};
                        }
                        bw.write(scores[3] + "");
                        bw.write(":");

                        // PL
                        bw.write(scores[0] + "," + scores[1] + "," + scores[2]);
                    } else {
                        bw.write(".:.,.");
                    }
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();;
            throw new IllegalArgumentException("Error writing VCF file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        }
        return filename;
    }

    /**
     * Writes given set of alignments to a set of Plink files
     *
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static String writeToPlink(GenotypeTable alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter MAPbw = null;
        BufferedWriter PEDbw = null;
        String mapFileName = Utils.addSuffixIfNeeded(filename, ".plk.map");
        String pedFileName = Utils.addSuffixIfNeeded(filename, ".plk.ped");
        try {
            MAPbw = new BufferedWriter(new FileWriter(mapFileName), 1000000);
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.chromosomeName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(alignment.siteName(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write("-9"); // genetic distance unavailable
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.chromosomalPosition(site))); // position
                MAPbw.write("\n");
            }
            MAPbw.close();

            PEDbw = new BufferedWriter(new FileWriter(pedFileName), 1000000);
            // Compiled : Pattern
            Pattern splitter = Pattern.compile(":");
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] name = splitter.split(alignment.taxaName(taxa).trim());
                if (name.length != 1) {
                    PEDbw.write(name[1]); // namelvl 1 if is available
                } else {
                    PEDbw.write("-9");
                }
                PEDbw.write(delimChar);
                PEDbw.write(alignment.taxaName(taxa).trim()); // namelvl 0
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
                    String[] b = getSNPValueForPlink(alignment.genotypeAsStringArray(taxa, site));
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
    public static String writeToFlapjack(GenotypeTable alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        String mapFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.map");
        String genoFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.geno");
        try {

            BufferedWriter MAPbw = new BufferedWriter(new FileWriter(mapFileName), 1000000);
            BufferedWriter DATbw = new BufferedWriter(new FileWriter(genoFileName), 1000000);
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.siteName(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write(alignment.chromosomeName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.chromosomalPosition(site))); // position
                MAPbw.write("\n");
                DATbw.write(delimChar);
                DATbw.write(alignment.siteName(site));
            }
            MAPbw.close();
            DATbw.write("\n");
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                DATbw.write(alignment.taxaName(taxa).trim());
                DATbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    String[] b = alignment.genotypeAsStringArray(taxa, site);
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

    public static String saveDelimitedAlignment(GenotypeTable theAlignment, String delimit, String saveFile) {

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
            int numSites = theAlignment.numberOfSites();
            for (int j = 0; j < numSites; j++) {
                bw.write(delimit);
                bw.write(String.valueOf(theAlignment.chromosomalPosition(j)));
            }
            bw.write("\n");

            for (int r = 0, n = theAlignment.numberOfTaxa(); r < n; r++) {
                bw.write(theAlignment.taxaName(r));
                for (int i = 0; i < numSites; i++) {
                    bw.write(delimit);
                    bw.write(theAlignment.genotypeAsString(r, i));
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
    public static void printSequential(GenotypeTable a, PrintWriter out) {
        // PHYLIP header line
        out.println("  " + a.numberOfTaxa() + " " + a.numberOfSites() + "  S");

        // Print sequences
        for (int s = 0; s < a.numberOfTaxa(); s++) {
            int n = 0;
            while (n < a.numberOfSites()) {
                if (n == 0) {
                    format.displayLabel(out, a.taxaName(s), 10);
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
    public static void printInterleaved(GenotypeTable a, PrintWriter out) {
        int n = 0;

        // PHYLIP header line
        out.println("  " + a.numberOfTaxa() + " " + a.numberOfSites());

        // Print sequences
        while (n < a.numberOfSites()) {
            for (int s = 0; s < a.numberOfTaxa(); s++) {
                if (n == 0) {
                    format.displayLabel(out, a.taxaName(s), 10);
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
    public static void printCLUSTALW(GenotypeTable a, PrintWriter out) {
        int n = 0;

        // CLUSTAL W header line
        out.println("CLUSTAL W multiple sequence alignment");
        out.println();

        // Print sequences
        while (n < a.numberOfSites()) {
            out.println();
            for (int s = 0; s < a.numberOfTaxa(); s++) {
                format.displayLabel(out, a.taxaName(s), 10);
                out.print("     ");

                printNextSites(a, out, false, s, n, 50);
                out.println();
            }
            // Blanks in status line are necessary for some parsers)
            out.println("               ");
            n += 50;
        }
    }

    private static void printNextSites(GenotypeTable a, PrintWriter out, boolean chunked, int seq, int start, int num) {
        // Print next num characters
        for (int i = 0; (i < num) && (start + i < a.numberOfSites()); i++) {
            // Chunks of 10 characters
            if (i % 10 == 0 && i != 0 && chunked) {
                out.print(' ');
            }
            out.print(a.genotypeAsString(seq, start + i));
        }
    }

    public static String writeAlignmentToSerialGZ(GenotypeTable sba, String outFile) {

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
