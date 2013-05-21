package net.maizegenetics.pd;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import net.maizegenetics.pal.alignment.BitNucleotideAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.util.Utils;

public class PDAnnotation {

    private static final String PHYSICAL_POSITIONS = HapMapHDF5Constants.POSITIONS;
    private static final String CHROMOSOME = "_C";
    private static final String GWAS_TRAIT = "_GT";
    private static final String MINOR_ALLELE_FREQUENCY = "MAF";
    private static final String MAJOR_ALLELE = "MajorAllele";
    private static final String MINOR_ALLELE = "MinorAllele";
    private static final String HAS_DATA = "HasData"; // summary index where if any trait has a value at that location, value is set to 1
    private boolean myIsSBit = true;
    private String[] chromosomes = new String[]{"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"};
    private String hapMapFile_prefix = "/maizeHapMapV2_B73RefGenV2_201203028_";
    private String hapMapFile_suffix = "h.hmp.txt.gz";
    private static final int PHYSICAL_POSITION_COLUMN = 0;
    private static final int MINOR_ALLELE_FREQUENCY_COLUMN = 1;
    private static final int COLUMN_OFFSET = 1; // one for physical position column and another for minor allele frequency column
    
    private int[][] allPositions;

    //private 
    public PDAnnotation(String hapMapPath, String pathGwas, String annoPath, String outputFile,
            int startChr, int endChr) {
        File aHapMapDir = new File(hapMapPath);
        File aGwasDir = new File(pathGwas);
        loadAllChromosomes(aHapMapDir, outputFile, startChr);
        loadGWAS(aGwasDir, outputFile, startChr);
    }



    public void loadAllChromosomes(File hapMapDir, String outputFile, int currChr) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(outputFile);
        config.overwrite();
        IHDF5Writer writer = config.writer();
        int[] hasData = null;  //recorder if there is any GWAS data for a site; TODO perhaps increment
        String chromosomeFile = hapMapDir + hapMapFile_prefix + "chr" + currChr + hapMapFile_suffix;
        System.out.println("Loading:" + chromosomeFile);
        Alignment bna = ImportUtils.readFromHapmap(chromosomeFile, myIsSBit, null /*progressListener*/);
        System.out.printf("Sites:%d StartPosition:%d EndPosition:%d %n", bna.getSiteCount(), bna.getPositionInLocus(0), bna.getPositionInLocus(bna.getSiteCount() - 1));
        bna.optimizeForSites(null);
        int siteCnt = bna.getSiteCount();
        int[] alignmentPhysPos = bna.getPhysicalPositions();
        hasData = new int[siteCnt];
        float[] maf = new float[siteCnt];
        byte[] mjAllele = new byte[siteCnt];
        byte[] mnAllele = new byte[siteCnt];
        for (int j = 0; j < siteCnt; j++) {
            mjAllele[j] = bna.getMajorAlleleAsString(j).getBytes()[0];
            mnAllele[j] = bna.getMinorAlleleAsString(j).getBytes()[0];
            maf[j] = (float) bna.getMinorAlleleFrequency(j);
        }
        //write positions to hdf "pos"+chromosome
        String chrGroup = "chr" + currChr + "/";
        writer.createGroup(chrGroup);
        writer.setIntAttribute(chrGroup, HapMapHDF5Constants.NUM_SITES, siteCnt);
        writer.createIntArray(chrGroup + HapMapHDF5Constants.POSITIONS, alignmentPhysPos.length);
        writer.writeIntArray(chrGroup + HapMapHDF5Constants.POSITIONS, alignmentPhysPos);

        //write alleles to hdf "allele"+chromosome
        // which version? String[][] ?
        writer.createByteArray(chrGroup + PDAnnotation.MAJOR_ALLELE, mjAllele.length);
        writer.writeByteArray(chrGroup + PDAnnotation.MAJOR_ALLELE, mjAllele);
        writer.createByteArray(chrGroup + PDAnnotation.MINOR_ALLELE, mnAllele.length);
        writer.writeByteArray(chrGroup + PDAnnotation.MINOR_ALLELE, mnAllele);
        // write minor allele frequencies
        writer.createFloatArray(chrGroup + PDAnnotation.MINOR_ALLELE_FREQUENCY, maf.length);
        writer.writeFloatArray(chrGroup + PDAnnotation.MINOR_ALLELE_FREQUENCY, maf);

        writer.createGroup(chrGroup + "GWAS");
        writer.createGroup(chrGroup + "GenomicAnno");
        writer.createGroup(chrGroup + "PopgenAnno");




        writer.close();
    }
 

    private void loadGWAS(File gwasFileIn, String outputFile, int currChr) {
        IHDF5Writer writer = HDF5Factory.open(outputFile);
//        FolderParser fp = new FolderParser(gwasDirIn);
        
        int traitIndex = 0; 
        int chrIndex = 1;
        int physPosIndex = 2;
        int resultIndex = 5;
        String[] traits =getGWASTraits(gwasFileIn, traitIndex, "\t");
        String chrGroup = "chr" + currChr + "/";
        //read in all chromosome position
        //create a method to hold this memory
        int[] positions = writer.readIntArray(chrGroup + HapMapHDF5Constants.POSITIONS);
        

        for (int j = 0; j < traits.length; j++) {
           
            System.out.println(gwasFileIn.toString());
            // pull out the physical location and p-value
            int posMatch = 0, posMisMatch = 0;
            float[] rmip = new float[positions.length];
            Arrays.fill(rmip, Float.NaN);
            try {
                BufferedReader fileIn = Utils.getBufferedReader(gwasFileIn, 1000000);
                String s;
                while ((s = fileIn.readLine()) != null) {
                    String[] fields = s.split("\t");
                    try {
                        int theChr = Integer.parseInt(fields[chrIndex]);
                        int position = Integer.parseInt(fields[physPosIndex]);
                        float rmipValue = Float.parseFloat(fields[resultIndex]);
                        if(theChr!=9) continue;
                        if(position>3600000) continue;
                        //int site = Arrays.binarySearch(allPositions[theChr-1], position);
                        int site = Arrays.binarySearch(positions, position); 
                        if (site < 0) {
                            System.out.println("Error Position not found:"+position);
                            posMisMatch++;
                        } else {
                            posMatch++;
                            rmip[site] = rmipValue;
                            System.out.printf("Hit Chr:%d Position:%d site:%d %n ",theChr,position, site);
                        }

                    } catch (Exception e) {
                        //                     System.out.println("Header");
                    }
                }
            } catch (IOException e) {
                System.out.println("IOError");
                e.printStackTrace();
            }
            System.out.printf("Position matches:%d errors:%d %n", posMatch, posMisMatch);
            String dataSetName = chrGroup + "GWAS/" + traits[j];
            writer.createFloatArray(dataSetName, rmip.length);
            writer.writeFloatArray(dataSetName, rmip);
        } // end of traits loop
    }
    
    private String[] getGWASTraits(File gwasResults, int traitIndex, String delimiter){
        BufferedReader br = Utils.getBufferedReader(gwasResults, 1000000);
        String line = null;
        Set<String> aSet = new HashSet();
        try{
            while((line =  br.readLine()) != null){
                String[] fields = line.split(delimiter);
                
                aSet.add(fields[traitIndex]);
            }
        }catch(IOException ioe){
            ioe.printStackTrace();
        }
    
        String[] result = new String[aSet.size()];
        aSet.toArray(result);
        return result;
    }

    public static void main(String[] args) {
//        String hapMapPath = "/local/workdir/dek29/pad/pad/HapMapV2RefGenV2/";
////        String pathGwas = "/local/workdir/dek29/pd/gwas_results/";
//        String PDfile = "/home/dek29/Documents/PolyDesc/output/testPD.h5";
//        String pathGwas = "/home/dek29/Documents/PolyDesc/20130521_fromJason/gwas_hits_all.txt";
        
         String hapMapPath = "/Volumes/LaCie/HapMapV2/compressed/";
        String pathGwas = "/Volumes/LaCie/PolymorphismDescriptors/gwas_hits_all.txt";
        String annoPath = "/Volumes/LaCie/HapMapV2/compressed/";
        String PDfile = "/Volumes/LaCie/PolymorphismDescriptors/testPD.h5";

        PDAnnotation p = new PDAnnotation(hapMapPath, pathGwas, annoPath, PDfile, 9, 9);
        //       p.init();
    }
}
