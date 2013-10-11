package net.maizegenetics.pal.alignment.io;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.ProjectionBuilder;
import net.maizegenetics.pal.alignment.genotype.ProjectionGenotype;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.pal.util.DonorHaplotypes;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.util.*;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class ProjectionAlignmentIO {

    public static Alignment getInstance(String paFile, String baseHighDensityAlignmentFile) {
        return getInstance(paFile, ImportUtils.readFromHapmap(baseHighDensityAlignmentFile, null));
    }

    public static Alignment getInstance(String paFile, Alignment baseHighDensityAlignment) {
        BufferedReader br = null;
        try {
            br = Utils.getBufferedReader(paFile);
            String[] sl=Utils.readLineSkipComments(br).split("\t");
            int baseTaxaCnt=Integer.parseInt(sl[0]);
            if(baseTaxaCnt!=baseHighDensityAlignment.getSequenceCount()) {
                System.err.println("Error in number of base taxa"); return null;
            }
            int taxaCnt=Integer.parseInt(sl[1]);
            for (int i = 0; i < baseTaxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                //change to hash map
                int index=Integer.parseInt(sl[0]);
                if(!baseHighDensityAlignment.getFullTaxaName(index).equals(sl[1])) {
                    System.err.println("Names or order does not agree with base taxa"); return null;
                }
            }
            SortedMap<Taxon,NavigableSet<DonorHaplotypes>> allBreakPoints=new TreeMap<>();
            for (int i = 0; i < taxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                int breakTotal=sl.length-1;
                if(breakTotal==0) continue;  //no data
                NavigableSet<DonorHaplotypes> breakForTaxon=new TreeSet<>();
                for (int bp = 0; bp < breakTotal; bp++) {
                    String[] bptext=sl[bp+1].split(":");
                    Chromosome chr=new Chromosome(bptext[0]);
                    DonorHaplotypes dh=new DonorHaplotypes(chr, Integer.parseInt(bptext[1]), Integer.parseInt(bptext[1]),
                            Integer.parseInt(bptext[2]), Integer.parseInt(bptext[3]));
                    breakForTaxon.add(dh);
                }
                allBreakPoints.put(new Taxon(sl[0]), breakForTaxon);
            }
            return ProjectionBuilder.getInstance(baseHighDensityAlignment, allBreakPoints);
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error reading Projection file: " + paFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void save(String outfile, Alignment pa) {
        if(!(pa.getGenotypeMatrix() instanceof ProjectionGenotype)) {
            throw new UnsupportedOperationException("Save only works for Alignments with projection genotypes");
        }
//        BufferedWriter bw = null;
//        try {
//            String fullFileName = Utils.addSuffixIfNeeded(outfile, ".pa.txt.gz", new String[]{".pa.txt", ".pa.txt.gz"});
//            bw = Utils.getBufferedWriter(fullFileName);
//            bw.write(myBaseAlignment.getSequenceCount()+"\t"+pa.getSequenceCount()+"\n");
//            bw.write("#Donor Haplotypes\n");
//            for (int i = 0; i < myBaseAlignment.getSequenceCount(); i++) {
//                bw.write(i+"\t"+myBaseAlignment.getFullTaxaName(i)+"\n");
//            }
//            bw.write("#Taxa Breakpoints\n");
//            bw.write("#Block are defined position:donor1:donor2 (-1 means no hypothesis)\n");
//            for (int i = 0; i < pa.getSequenceCount(); i++) {
//                bw.write(pa.getFullTaxaName(i)+"\t");
////                for (int p = 0; (myPosBreaks[i]!=null)&&(p < myPosBreaks[i].length); p++) {
////                    bw.write(myPosBreaks[i][p]+":"+myHDTaxa[i][p][0]+":"+myHDTaxa[i][p][1]+"\t");
////                }
//                bw.write("\n");
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//            throw new IllegalArgumentException("Error writing Projection file: " + outfile + ": " + ExceptionUtils.getExceptionCauses(e));
//        } finally {
//            try {bw.close();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }
    }

}
