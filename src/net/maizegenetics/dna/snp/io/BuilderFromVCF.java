package net.maizegenetics.dna.snp.io;

import com.google.common.base.Splitter;
import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.depth.AlleleDepthBuilder;
import net.maizegenetics.dna.snp.depth.AlleleDepthUtil;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
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
    private static final Pattern WHITESPACE_PATTERN=Pattern.compile("[\\s]+");
    private HeaderPositions hp=null;
    private final String infile;
    private boolean includeDepth=false;

    private BuilderFromVCF(String infile) {
        this.infile=infile;
    }

    public static BuilderFromVCF getBuilder(String infile) {
        return new BuilderFromVCF(infile);
    }

    public BuilderFromVCF keepDepth() {
        includeDepth=true;
        return this;
    }

    //TODO provide options on caching to use, read only some sites, etc.
    public GenotypeTable build() {
        long time=System.nanoTime();
        GenotypeTable result=null;
        try {
            int numThreads=Runtime.getRuntime().availableProcessors();
            ExecutorService pool=Executors.newFixedThreadPool(numThreads);
            BufferedReader r=Utils.getBufferedReader(infile, -1);
            //Read the ## annotation rows
            String currLine;
            Map<String,String> infoMap=new HashMap<>();
            Map<String,String> formatMap=new HashMap<>();
            Map<String,SetMultimap<String,String>> sampAnnoBuild=new TreeMap<>();
            while (((currLine=r.readLine())!=null)&&(currLine.startsWith("##"))) {
                String[] cat=currLine.split("=",2);
                if(cat.length<2) continue;
                switch (cat[0]) {
//                    case "##INFO":
//                        infoMap.put(mapOfAnno.get("ID"), mapOfAnno.get("Description"));
//                        break;
//                    case "##FILTER":break;
//                    case "##FORMAT":
//                        formatMap.put(mapOfAnno.get("ID"),mapOfAnno.get("Description"));
//                        break;
                    case "##SAMPLE":
                        SetMultimap<String, String> mapOfAnno=TaxaListIOUtils.parseVCFHeadersIntoMap(cat[1]);
                        String taxaID=mapOfAnno.get("ID").iterator().next();
                        if(taxaID==null) break;
                        sampAnnoBuild.put(taxaID,mapOfAnno);
                        break;
                    case "##PEDIGREE":break;
                    default : break;
                }

//                System.out.println(currLine);
            }
            TaxaList taxaList=processTaxa(currLine,sampAnnoBuild);
            int linesAtTime=1<<12;  //this is a critical lines with 20% or more swings.  Needs to be optimized with transposing
            //  int linesAtTime=1<<8;  //better for with lots of taxa.
            ArrayList<String> txtLines=new ArrayList<>(linesAtTime);
            ArrayList<ProcessVCFBlock> pbs=new ArrayList<>();
            int lines=0;
            while ((currLine=r.readLine())!=null) {
                if(currLine.startsWith("#")) continue;
                txtLines.add(currLine);
                lines++;
                if (lines%linesAtTime==0) {
                    ProcessVCFBlock pb=ProcessVCFBlock.getInstance(taxaList.numberOfTaxa(), hp, txtLines);
                    pbs.add(pb);
                         pb.run();
                    //pool.execute(pb);
                    txtLines=new ArrayList<>(linesAtTime);
                }
            }
            r.close();
            if (txtLines.size()>0) {
                ProcessVCFBlock pb=ProcessVCFBlock.getInstance(taxaList.numberOfTaxa(), hp, txtLines);
                pbs.add(pb);
                pool.execute(pb);
            }
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BuilderFromHapMap: processing threads timed out.");
            }
            int currentSite=0;
            PositionListBuilder posBuild=new PositionListBuilder();
            GenotypeCallTableBuilder gb=GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.numberOfTaxa(), lines);
            AlleleDepthBuilder db=null;
            if(includeDepth) db=AlleleDepthBuilder.getInstance(taxaList.numberOfTaxa(),lines,6);
            for (ProcessVCFBlock pb : pbs) {
                posBuild.addAll(pb.getBlkPosList());
                byte[][] bgTS=pb.getGenoTS();
                for (int t=0; t<bgTS.length; t++) {
                    gb.setBaseRangeForTaxon(t, currentSite, bgTS[t]);
                }
                if(includeDepth) {
                    byte[][][] bdTS=pb.getDepthTS();
                    for (int t=0; t<bgTS.length; t++) {
                        db.setDepthRangeForTaxon(t, currentSite, bdTS[t]);
                    }
                }
                currentSite+=pb.getSiteNumber();
            }
            if (posBuild.validateOrdering()==false) {
                throw new IllegalStateException("BuilderFromHapMap: Ordering incorrect HapMap must be ordered by position");
            }
            GenotypeCallTable g=gb.build();
            if(includeDepth) {result=GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList, null, db.build());}
            else {result=GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList);}
        } catch (IOException|InterruptedException e) {
            e.printStackTrace();
        }
        long totalTime=System.nanoTime()-time;
        System.out.printf("BuilderFromVCF data timing %gs %n", totalTime/1e9);
        return result;
    }

//    static SetMultimap<String,String> parseVCFHeadersIntoMap(String s) {
//        if(s==null) return null;
//        if(!(s.startsWith("<") && s.endsWith(">"))) return null;
//        String value=s.substring(1,s.length()-1);
// //       System.out.println(s);
//        value=getReplaceCommaWithinQuote(value);
//
//        ImmutableSetMultimap.Builder<String,String> im=new ImmutableSetMultimap.Builder<String,String>()
//                .orderKeysBy(Ordering.natural()).orderValuesBy(Ordering.natural());
//        for (String s1 : Splitter.on(",").trimResults().split(value)) {
//            String[] ssEntry=s1.split("=",2);
//            String v=ssEntry[1];
//            if(v.contains(""+(char) ((int) ','+256)) || v.contains(""+(char) ((int) '='+256))) {
//                v=v.replace((char)((int)','+256),',');
//                v=v.replace((char)((int)'='+256),'=');
//            }
//            im.put(ssEntry[0],v);
//        }
//         return im.build();
//    }

    private static String getReplaceCommaWithinQuote(String s) {
        StringBuilder sb =new StringBuilder(s);
        boolean inQuote=false;
        for (int i=0; i<sb.length(); i++) {
            if(sb.charAt(i)=='\"') inQuote=(!inQuote);
            if(inQuote && sb.charAt(i)==',') sb.setCharAt(i,(char)((int)','+256));//(char)167);
            if(inQuote && sb.charAt(i)=='=') sb.setCharAt(i,(char)((int)'='+256));//(char)167);
        }
        return sb.toString();
    }

    private TaxaList processTaxa(String readLn, Map<String,SetMultimap<String,String>> taxaAnnotation) {
        String[] header=WHITESPACE_PATTERN.split(readLn);
        hp=new HeaderPositions(header);
        int numTaxa=header.length-hp.NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i=0; i<numTaxa; i++) {
            String taxonID=header[i+hp.NUM_HAPMAP_NON_TAXA_HEADERS];
            Taxon.Builder at=new Taxon.Builder(taxonID);
            SetMultimap<String,String> taMap=taxaAnnotation.get(taxonID);
            if(taMap!=null) {
                for (Map.Entry<String,String> en : taMap.entries()) {
                    if(en.getKey().equals("ID")) continue; //skip the IDs as these became the name
                    String s=en.getValue().replace("\"","");
                    at.addAnno(en.getKey(),s);
                }
            }
            tlb.add(at.build());
        }
        return tlb.build();
    }




}

class HeaderPositions {
    final int NUM_HAPMAP_NON_TAXA_HEADERS;
    final int GENOIDX;
    final int SNPID_INDEX;
  //  final int VARIANT_INDEX;
    final int FILTER_INDEX;
    final int QUAL_INDEX;
    final int CHROMOSOME_INDEX;
    final int POSITION_INDEX;
    final int REF_INDEX;
    final int ALT_INDEX;
    final int INFO_INDEX;
    final int FORMAT_INDEX;

    public HeaderPositions(String[] header){
        int chrIdx=firstEqualIndex(header,"#CHROM");
        if(chrIdx<0) chrIdx=firstEqualIndex(header,"#CHR");
        CHROMOSOME_INDEX=chrIdx;
        POSITION_INDEX=firstEqualIndex(header,"POS");
        SNPID_INDEX=firstEqualIndex(header,"ID");
        REF_INDEX=firstEqualIndex(header,"REF");
        ALT_INDEX=firstEqualIndex(header,"ALT");
        QUAL_INDEX=firstEqualIndex(header,"QUAL");
        FILTER_INDEX=firstEqualIndex(header,"FILTER");
        INFO_INDEX=firstEqualIndex(header,"INFO");
        FORMAT_INDEX=firstEqualIndex(header,"FORMAT");

        NUM_HAPMAP_NON_TAXA_HEADERS=Math.max(INFO_INDEX,FORMAT_INDEX)+1;
        GENOIDX=NUM_HAPMAP_NON_TAXA_HEADERS;
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i=0; i<sa.length; i++) {
            if(sa[i].equals(match)) return i;
        }
        return -1;
    }

}

class ProcessVCFBlock implements Runnable {

    private static final Pattern WHITESPACE_PATTERN=Pattern.compile("\\s");
    private static final Pattern SLASH_PATTERN=Pattern.compile("/");
//    private static final int NUM_HAPMAP_NON_TAXA_HEADERS=5;
//    private static final int GENOIDX=NUM_HAPMAP_NON_TAXA_HEADERS;
//    private static final int SNPID_INDEX=-1;
//    private static final int VARIANT_INDEX=-1;
//    private static final int CHROMOSOME_INDEX=0;
//    private static final int POSITION_INDEX=1;
//    private static final int REF_INDEX=2;
//    private static final int ALT_INDEX=3;
//    private static final int INFO_INDEX=3;
    private final HeaderPositions hp;
    private final int taxaN;
    private final int siteN;
    private ArrayList<String> txtL;
    private byte[][] gTS;  //genotypes
    private byte[][][] dTS; //depth
    private final ArrayList<Position> blkPosList;


    private ProcessVCFBlock(int taxaN, HeaderPositions hp, ArrayList<String> txtL) {
        this.taxaN=taxaN;
        this.siteN=txtL.size();
        this.txtL=txtL;
        this.hp=hp;
        blkPosList=new ArrayList<>(siteN);

    }

    public static ProcessVCFBlock getInstance(int taxaN, HeaderPositions hp, ArrayList<String> txtL) {
        return new ProcessVCFBlock(taxaN, hp, txtL);
    }

    //@Override
    public void run() {
        Map<String, Chromosome> chromosomeLookup=new HashMap<>();
        gTS=new byte[taxaN][siteN];
        dTS=new byte[taxaN][6][siteN];
        for (int s=0; s<siteN; s++) {
            //really needs to use a Splitter iterator to make this cleaner if it is performant
            String input=txtL.get(s);
            try{
                int[] tabPos=new int[hp.NUM_HAPMAP_NON_TAXA_HEADERS+taxaN];
                int tabIndex=0;
                int len=input.length();
                for (int i=0; (tabIndex<hp.NUM_HAPMAP_NON_TAXA_HEADERS+taxaN)&&(i<len); i++) {
                    if (input.charAt(i)=='\t') {
                        tabPos[tabIndex++]=i;
                    }
                }
                String chrName=input.substring(0, tabPos[hp.CHROMOSOME_INDEX]);
                Chromosome currChr=chromosomeLookup.get(chrName);
                if (currChr==null) {
                    currChr=new Chromosome(new String(chrName));
                    chromosomeLookup.put(chrName, currChr);
                }
                String snpID=null;
                if(hp.SNPID_INDEX>0) snpID=input.substring(tabPos[hp.SNPID_INDEX-1]+1, tabPos[hp.SNPID_INDEX]);
                String refS=input.substring(tabPos[hp.REF_INDEX-1]+1, tabPos[hp.REF_INDEX]);
                String alt=input.substring(tabPos[hp.ALT_INDEX-1]+1, tabPos[hp.ALT_INDEX]);
                String variants;
                if(alt.equals(".")) {variants=refS;}
                else {variants=(refS+"/"+alt).replace(',','/')
                        .replace("<INS>", "+").replace('I', '+')
                        .replace("<DEL>", "-").replace('D', '-');}
                GeneralPosition.Builder apb=new GeneralPosition.Builder(currChr, Integer.parseInt(input.substring(tabPos[hp.POSITION_INDEX-1]+1, tabPos[hp.POSITION_INDEX])))
                        .knownVariants(variants) //TODO strand, variants,
                        ;
                if(snpID!=null && !snpID.equals(".")) {
                    apb.snpName(snpID);
                }
                byte[] alleles=new byte[(variants.length()+1)/2];
                for (int i = 0, varInd=0; i < alleles.length; i++, varInd+=2) {
                    alleles[i]=NucleotideAlignmentConstants.getNucleotideAlleleByte(variants.charAt(varInd));
                }
                apb.allele(Position.Allele.REF, alleles[0]);
                for(String annoS: Splitter.on(";").split(input.substring(tabPos[hp.INFO_INDEX-1]+1, tabPos[hp.INFO_INDEX]))) {
                    apb.addAnno(annoS);
                }
                blkPosList.add(apb.build());
                final int iGT=0; //genotype index
                int iAD=-1,iDP=-1,iGQ=-1, iPL=-1;  //alleleDepth, overall depth, genotypeQuality, phredGenotypeLikelihoods
                if(hp.FORMAT_INDEX>=0) {
                    String[] formatS=input.substring(tabPos[hp.FORMAT_INDEX-1]+1, tabPos[hp.FORMAT_INDEX]).split(":");
                    iAD=firstEqualIndex(formatS,"AD");
                }
                int t=0;
                for(String taxaAllG: Splitter.on("\t").split(input.substring(tabPos[hp.NUM_HAPMAP_NON_TAXA_HEADERS-1]+1))) {
                    int f=0;
                    for(String fieldS: Splitter.on(":").split(taxaAllG)) {
                        if(f==iGT) {
                            int a1=fieldS.charAt(0)-'0';
                            int a2=fieldS.charAt(2)-'0';
                            if(a1<0 || a2<0 ) {gTS[t][s]=GenotypeTable.UNKNOWN_DIPLOID_ALLELE;}
                            else {gTS[t][s]=GenotypeTableUtils.getDiploidValue(alleles[a1],alleles[a2]);}
                        } else if(f==iAD) {
                            int i=0;
                            for(String ad: Splitter.on(",").split(fieldS)){
                                if(alleles[i]==GenotypeTable.UNKNOWN_ALLELE) {  //no position for depth of unknown alleles, so skip
                                    i++;
                                    continue;}
                                int adInt=Integer.parseInt(ad);
                                dTS[t][alleles[i++]][s]=AlleleDepthUtil.depthIntToByte(adInt);
                            }
                        }
                        f++;
                    }
                    t++;
                }
            } catch(Exception e) {
                System.err.println("Err Site Number:"+s);
                System.err.println("Err:"+input);
                throw e;
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

    byte[][][] getDepthTS() {
        return dTS;
    }

    ArrayList<Position> getBlkPosList() {
        return blkPosList;
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i=0; i<sa.length; i++) {
            if(sa[i].equals(match)) return i;
        }
        return -1;
    }
}

