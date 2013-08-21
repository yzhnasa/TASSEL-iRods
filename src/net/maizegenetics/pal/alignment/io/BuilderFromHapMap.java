package net.maizegenetics.pal.alignment.io;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.site.*;
import net.maizegenetics.pal.taxa.AnnotatedTaxon;
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
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class BuilderFromHapMap {
    private static final Logger myLogger = Logger.getLogger(BuilderFromHapMap.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private final String infile;

    private BuilderFromHapMap(String infile) {
        this.infile=infile;
    }

    public static BuilderFromHapMap getBuilder(String infile) {
        return new BuilderFromHapMap(infile);
    }

    //TODO provide options on caching to use, read only some sites, etc.

    public AlignmentNew build() {
        long time=System.nanoTime();
        AlignmentNew result=null;
        try {
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);
            BufferedReader r=Utils.getBufferedReader(infile, -1);
            TaxaList taxaList=processTaxa(r.readLine());
            String currLine;
            int linesAtTime=1<<12;
            ArrayList<String> txtLines=new ArrayList<>(linesAtTime);
            ArrayList<ProcessBlock> pbs=new ArrayList<>();
            int lines=0;
            while((currLine=r.readLine())!=null) {
                txtLines.add(currLine);
                lines++;
                if(lines%linesAtTime==0) {
                    ProcessBlock pb=ProcessBlock.getInstance(pbs.size(),taxaList.getTaxaCount(),txtLines);
                    pbs.add(pb);
               //     pb.run();
                    pool.execute(pb);
                    txtLines=new ArrayList<>(linesAtTime);
                }
            }
            r.close();
            if(txtLines.size()>0) {
                ProcessBlock pb=ProcessBlock.getInstance(pbs.size(),taxaList.getTaxaCount(),txtLines);
                pbs.add(pb);
                pool.execute(pb);
            }
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
            }
            int currentSite=0;
            AnnotatedPositionArrayList.Builder posBuild=new AnnotatedPositionArrayList.Builder();
            GenotypeBuilder gb = GenotypeBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.getTaxaCount(), lines);
            for(ProcessBlock pb: pbs) {
                posBuild.addAll(pb.getBlkPosList());
                byte[][] bgTS=pb.getGenoTS();
                for (int t=0; t<bgTS.length; t++) {
                    gb.setBaseRangeForTaxon(t,currentSite,bgTS[t]);
                }
                currentSite+=pb.getSiteNumber();
            }
            Genotype g = gb.build();
            result=new CoreAlignment(g,posBuild.build(),taxaList,null,null);
        } catch (IOException|InterruptedException e) {
            e.printStackTrace();
        }
        long totalTime=System.nanoTime()-time;
        System.out.printf("ImportUtil ReadText data timing %gs %n", totalTime/1e9);
        return result;
    }

    private static TaxaList processTaxa(String readLn) {
        String[] header = WHITESPACE_PATTERN.split(readLn);
        int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (int i = 0; i < numTaxa; i++) {
            AnnotatedTaxon at = new AnnotatedTaxon.Builder(header[i + NUM_HAPMAP_NON_TAXA_HEADERS])
                    .build();
            tlb.add(at);
        }
        return tlb.build();
    }

}

class ProcessBlock implements Runnable {
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private static final int HAPMAP_SNPID_COLUMN_INDEX = 0;
    private static final int HAPMAP_CHROMOSOME_COLUMN_INDEX = 2;
    private static final int HAPMAP_POSITION_COLUMN_INDEX = 3;
    private final int order;
    private final int taxaN;
    private final int siteN;
    private ArrayList<String> txtL;
    private byte[][] gTS;
    private final ArrayList<AnnotatedPosition> blkPosList;
    private final byte[] convert;

    private ProcessBlock(int order, int taxaN, ArrayList<String> txtL) {
        this.order=order;
        this.taxaN=taxaN;
        this.siteN=txtL.size();
        this.txtL=txtL;
        blkPosList=new ArrayList<>(siteN);

        convert = new byte[128];
        for (int i = 0; i < convert.length; i++) {
            try {
                convert[i] = NucleotideAlignmentConstants.getNucleotideDiploidByte((char) i);
            } catch (IllegalArgumentException e) {
                convert[i] = Alignment.UNKNOWN_DIPLOID_ALLELE;
            }
        }
    }

    public static ProcessBlock getInstance(int order, int taxaN, ArrayList<String> txtL) {
        return new ProcessBlock(order,taxaN, txtL);
    }

    @Override
    public void run() {
        Map<String, Chromosome> chromosomeLookup = new HashMap<>();
        gTS=new byte[taxaN][siteN];
        for (int s=0; s<siteN; s++) {
            String input=txtL.get(s);
            String[] tokens = WHITESPACE_PATTERN.split(input,NUM_HAPMAP_NON_TAXA_HEADERS+1);
//            int sumLength = NUM_HAPMAP_NON_TAXA_HEADERS;
//            for (int i = 0; i < NUM_HAPMAP_NON_TAXA_HEADERS; i++) {
//                sumLength += tokens[i].length();
//            }
            Chromosome currChr = chromosomeLookup.get(tokens[HAPMAP_CHROMOSOME_COLUMN_INDEX]);
            if (currChr == null) {
                currChr = new Chromosome(tokens[HAPMAP_CHROMOSOME_COLUMN_INDEX]);
                chromosomeLookup.put(tokens[HAPMAP_CHROMOSOME_COLUMN_INDEX], currChr);
            }
            CorePosition cp = new CorePosition.Builder(currChr, Integer.parseInt(tokens[HAPMAP_POSITION_COLUMN_INDEX]))
                    .snpName(tokens[HAPMAP_SNPID_COLUMN_INDEX])
                            //TODO                    strand, variants,
                    .build();
            CoreAnnotatedPosition apb = new CoreAnnotatedPosition.Builder(cp).build();
            blkPosList.add(apb);
            if(s==0) System.out.println(apb.getChromosome()+":"+tokens[11].length()/2);
            int avg=Math.round((float)tokens[11].length()/(float)taxaN);
            //System.out.println(tokens[11].split("\t").length);
            if(avg==2) {
                for (int i = 0; i < tokens[11].length(); i+=2) {
                    //System.out.println(i+":"+tokens[11].charAt(i));
                    gTS[i/2][s] = convert[tokens[11].charAt(i)];
                }
            } else if(avg==3) {
                for (int i = 0; i < tokens[11].length(); i+=3) {
                    //System.out.println(i+":"+tokens[11].charAt(i)+tokens[11].charAt(i));
                    gTS[i/3][s] =AlignmentUtils.getDiploidValue(convert[tokens[11].charAt(i)], convert[tokens[11].charAt(i+1)]);
                }
            }
        }
        txtL=null;
    }

    int getSiteNumber() {
        return siteN;
    }

    byte[][] getGenoTS() {
        return gTS;
    }

    ArrayList<AnnotatedPosition> getBlkPosList() {
        return blkPosList;
    }
}
