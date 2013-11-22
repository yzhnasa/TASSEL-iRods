package net.maizegenetics.dna.snp.io;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.AlignmentBuilder;
import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.dna.snp.AlignmentUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotype.Genotype;
import net.maizegenetics.dna.snp.genotype.GenotypeBuilder;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

/**
 * Create an alignment based on VCF format file (either .txt or compressed).  Alleles are set as global reference.
 * e.g. code <p></p>
 * {@code
 * Alignment a=BuilderFromVCF.getBuilder(infileName).build();
 * }
 * <p></p>
 * TODO:  Add filtering while reading, provide an option to define the alleles as reference and alternate
 *
 * @author Ed Buckler
 */
public class BuilderFromVCF {

    private static final Logger myLogger=Logger.getLogger(BuilderFromHapMap.class);
    private static final Pattern WHITESPACE_PATTERN=Pattern.compile("\\s");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS=5;
    private final String infile;

    private BuilderFromVCF(String infile) {
        this.infile=infile;
    }

    public static BuilderFromVCF getBuilder(String infile) {
        return new BuilderFromVCF(infile);
    }

    //TODO provide options on caching to use, read only some sites, etc.
    public Alignment build() {
        long time=System.nanoTime();
        Alignment result=null;
        try {
            int numThreads=Runtime.getRuntime().availableProcessors();
            ExecutorService pool=Executors.newFixedThreadPool(numThreads);
            BufferedReader r=Utils.getBufferedReader(infile, -1);
            TaxaList taxaList=processTaxa(r.readLine());
            String currLine;
            int linesAtTime=1<<12;  //this is a critical lines with 20% or more swings.  Needs to be optimized with transposing
            //  int linesAtTime=1<<8;  //better for with lots of taxa.
            ArrayList<String> txtLines=new ArrayList<>(linesAtTime);
            ArrayList<ProcessVCFBlock> pbs=new ArrayList<>();
            int lines=0;
            while ((currLine=r.readLine())!=null) {
                txtLines.add(currLine);
                lines++;
                if (lines%linesAtTime==0) {
                    ProcessVCFBlock pb=ProcessVCFBlock.getInstance(taxaList.getTaxaCount(), txtLines);
                    pbs.add(pb);
                    //     pb.run();
                    pool.execute(pb);
                    txtLines=new ArrayList<>(linesAtTime);
                }
            }
            r.close();
            if (txtLines.size()>0) {
                ProcessVCFBlock pb=ProcessVCFBlock.getInstance(taxaList.getTaxaCount(), txtLines);
                pbs.add(pb);
                pool.execute(pb);
            }
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BuilderFromHapMap: processing threads timed out.");
            }
            int currentSite=0;
            PositionListBuilder posBuild=new PositionListBuilder();
            GenotypeBuilder gb=GenotypeBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.getTaxaCount(), lines);
            for (ProcessVCFBlock pb : pbs) {
                posBuild.addAll(pb.getBlkPosList());
                byte[][] bgTS=pb.getGenoTS();
                for (int t=0; t<bgTS.length; t++) {
                    gb.setBaseRangeForTaxon(t, currentSite, bgTS[t]);
                }
                currentSite+=pb.getSiteNumber();
            }
            if (posBuild.validateOrdering()==false) {
                throw new IllegalStateException("BuilderFromHapMap: Ordering incorrect HapMap must be ordered by position");
            }
            Genotype g=gb.build();
            result=AlignmentBuilder.getInstance(g, posBuild.build(), taxaList);
        } catch (IOException|InterruptedException e) {
            e.printStackTrace();
        }
        long totalTime=System.nanoTime()-time;
        System.out.printf("BuilderFromHapMap data timing %gs %n", totalTime/1e9);
        return result;
    }

    private static TaxaList processTaxa(String readLn) {
        String[] header=WHITESPACE_PATTERN.split(readLn);
        int numTaxa=header.length-NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i=0; i<numTaxa; i++) {
            Taxon at=new Taxon.Builder(header[i+NUM_HAPMAP_NON_TAXA_HEADERS]).build();
            tlb.add(at);
        }
        return tlb.build();
    }
}

class ProcessVCFBlock implements Runnable {

    private static final Pattern WHITESPACE_PATTERN=Pattern.compile("\\s");
    private static final Pattern SLASH_PATTERN=Pattern.compile("/");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS=5;
    private static final int GENOIDX=NUM_HAPMAP_NON_TAXA_HEADERS;
    private static final int SNPID_INDEX=-1;
    private static final int VARIANT_INDEX=-1;
    private static final int CHROMOSOME_INDEX=0;
    private static final int POSITION_INDEX=1;
    private static final int REF_INDEX=2;
    private static final int ALT_INDEX=3;
    private static final int INFO_INDEX=3;
    private final int taxaN;
    private final int siteN;
    private ArrayList<String> txtL;
    private byte[][] gTS;
    private final ArrayList<Position> blkPosList;


    private ProcessVCFBlock(int taxaN, ArrayList<String> txtL) {
        this.taxaN=taxaN;
        this.siteN=txtL.size();
        this.txtL=txtL;
        blkPosList=new ArrayList<>(siteN);

    }

    public static ProcessVCFBlock getInstance(int taxaN, ArrayList<String> txtL) {
        return new ProcessVCFBlock(taxaN, txtL);
    }

    //@Override
    public void run() {
        Map<String, Chromosome> chromosomeLookup=new HashMap<>();
        gTS=new byte[taxaN][siteN];
        for (int s=0; s<siteN; s++) {
            String input=txtL.get(s);
            int[] tabPos=new int[NUM_HAPMAP_NON_TAXA_HEADERS+taxaN];
            int tabIndex=0;
            int len=input.length();
            for (int i=0; (tabIndex<NUM_HAPMAP_NON_TAXA_HEADERS+taxaN)&&(i<len); i++) {
                if (input.charAt(i)=='\t') {
                    tabPos[tabIndex++]=i;
                }
            }
            String chrName=input.substring(0, tabPos[CHROMOSOME_INDEX]);
            Chromosome currChr=chromosomeLookup.get(chrName);
            if (currChr==null) {
                currChr=new Chromosome(new String(chrName));
                chromosomeLookup.put(chrName, currChr);
            }
            String refS=input.substring(tabPos[REF_INDEX-1]+1, tabPos[REF_INDEX]);
            String alt=input.substring(tabPos[ALT_INDEX-1]+1, tabPos[ALT_INDEX]);
            String variants=(refS+"/"+alt).replace(',','/').replace('I','+').replace('D','-');
            GeneralPosition.Builder apb=new GeneralPosition.Builder(currChr, Integer.parseInt(input.substring(tabPos[POSITION_INDEX-1]+1, tabPos[POSITION_INDEX])))
                    .knownVariants(variants) //TODO                    strand, variants,
                    ;
            byte[] alleles=new byte[(variants.length()+1)/2];
            for (int i = 0, varInd=0; i < alleles.length; i++, varInd+=2) {
                alleles[i]=NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(varInd));
            }
            apb.allele(Position.Allele.REF, alleles[0]);
            blkPosList.add(apb.build());
//            System.out.println(s);
//            System.out.println(input);
            for (int t = 0; t < taxaN; t++) {
                int offset=tabPos[NUM_HAPMAP_NON_TAXA_HEADERS-1+t]+1;
                int a1=input.charAt(offset)-'0';
                int a2=input.charAt(offset+2)-'0';
                if(a1<0 || a2<0 ) {gTS[t][s]=Alignment.UNKNOWN_DIPLOID_ALLELE;}
                else {gTS[t][s]=AlignmentUtils.getDiploidValue(alleles[a1],alleles[a2]);}
              //  System.out.print(gTS[t][s] + "\t");
            }
//            System.out.println();
        }
        txtL=null;
    }

    int getSiteNumber() {
        return siteN;
    }

    byte[][] getGenoTS() {
        return gTS;
    }

    ArrayList<Position> getBlkPosList() {
        return blkPosList;
    }
}

