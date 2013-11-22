package net.maizegenetics.dna.snp.io;

import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.dna.snp.AlignmentBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.ProjectionBuilder;
import net.maizegenetics.dna.snp.genotype.ProjectionGenotype;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.util.DonorHaplotypes;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

/**
 * Methods for reading and writing ProjectionAlignments to files.
 *
 * @author Ed Buckler
 */
public class ProjectionAlignmentIO {

    public static Alignment getInstance(String paFile, String baseHighDensityAlignmentFile) {
        if(baseHighDensityAlignmentFile.endsWith(".h5")) return getInstance(paFile, AlignmentBuilder.getInstance(baseHighDensityAlignmentFile));
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
            Map<Integer,Integer> paIndexToBaseIndex=new HashMap<>();
            for (int i = 0; i < baseTaxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                int index=Integer.parseInt(sl[0]);
                Taxon taxon=new Taxon(sl[1]);
                List<Integer> matches=baseHighDensityAlignment.getTaxaList().getIndicesMatchingTaxon(taxon);
                if(matches.size()==0) {
                    throw new NoSuchElementException("Taxon "+sl[1]+" not found within base taxa");
                }
                if(matches.size()>1) {
                    throw new NoSuchElementException("Taxon "+sl[1]+" found multiple times within base taxa");
                }
                paIndexToBaseIndex.put(index, matches.get(0));
            }
            ProjectionBuilder pb=new ProjectionBuilder(baseHighDensityAlignment);
            for (int i = 0; i < taxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                int breakTotal=sl.length-1;
                if(breakTotal==0) continue;  //no data
                NavigableSet<DonorHaplotypes> breakForTaxon=new TreeSet<>();
                for (int bp = 0; bp < breakTotal; bp++) {
                    String[] bptext=sl[bp+1].split(":");
                    Chromosome chr=new Chromosome(bptext[0]);
                    int baseParent1=paIndexToBaseIndex.get(Integer.parseInt(bptext[3]));
                    int baseParent2=paIndexToBaseIndex.get(Integer.parseInt(bptext[4]));
                    DonorHaplotypes dh=new DonorHaplotypes(chr, Integer.parseInt(bptext[1]), Integer.parseInt(bptext[2]),
                            baseParent1, baseParent2);
                    breakForTaxon.add(dh);
                }
                pb.addTaxon(new Taxon(sl[0]), breakForTaxon);
            }
            return pb.build();
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

    public static void writeToFile(String outfile, Alignment pa) {
        if(!(pa.getGenotypeMatrix() instanceof ProjectionGenotype)) {
            throw new UnsupportedOperationException("Save only works for Alignments with projection genotypes");
        }
        ProjectionGenotype pg=(ProjectionGenotype)pa.getGenotypeMatrix();
        Alignment baseAlignment=pg.getBaseAlignment();
        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(outfile, ".pa.txt.gz", new String[]{".pa.txt", ".pa.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
            bw.write(baseAlignment.getSequenceCount()+"\t"+pa.getSequenceCount()+"\n");
            bw.write("#Donor Haplotypes\n");
            for (int i = 0; i < baseAlignment.getSequenceCount(); i++) {
                bw.write(i+"\t"+baseAlignment.getTaxaName(i)+"\n");
            }
            bw.write("#Taxa Breakpoints\n");
            bw.write("#Block are defined chr:startPos:endPos:donor1:donor2 (-1 means no hypothesis)\n");
            for (int i = 0; i < pa.getSequenceCount(); i++) {
                bw.write(pa.getTaxaName(i)+"\t");
                NavigableSet<DonorHaplotypes> theDH=pg.getDonorHaplotypes(i);
                for (DonorHaplotypes dh : theDH) {
                    bw.write(dh.getChromosome().getName()+":"+dh.getStartPosition()+":"+dh.getEndPosition()+":"+
                            dh.getParent1index()+":"+dh.getParent2index()+"\t");
                }
                bw.write("\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing Projection file: " + outfile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

}
