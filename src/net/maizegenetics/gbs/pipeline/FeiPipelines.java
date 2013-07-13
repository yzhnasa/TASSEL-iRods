/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.TagsOnGeneticMap;
import net.maizegenetics.gbs.maps.TagsOnPhysMapHDF5;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;

/**
 *
 * @author Fei Lu
 */
public class FeiPipelines {
    
    public FeiPipelines () {
        this.pipelinePE();
        //this.testPipeline();
    }
    
    public static void main (String[] args) {
        new FeiPipelines();
        
    }
    
    public void testPipeline () {
        //this.mkSmallTagCount();
        //this.mkFastq();
        //this.mkFasta();
        //this.mkTOPM();
        //this.mkTOPMHDF5();
    }
    
    public void mkTOPMHDF5 () {
        String infileS = "M:/GBStest/topm/bowtie2.topm.bin";
        String outfileS = "M:/GBStest/topm/bowtie2.topm.h5";
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap (infileS, false);
        TagsOnPhysMapHDF5.createFile(topm, outfileS, 4, 16);
    }
    
    public void mkTOPM () {
        String inputFileS = "M:/GBStest/alignment/bowtie2.sam";
        String outputFileS = "M:/GBStest/topm/bowtie2.topm.bin";
        new FeiUtils ().convertSAM2TOPM(inputFileS, outputFileS);
    }
    
    /**
     * Make FASTA file and do alignment using Blast, -m 8 -e 1e-10
     */
    public void mkFasta () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/alignment/small.fa";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        tc.toFASTA(outputFileS);
    }
    
    /**
     * Make Fastq file and do alignment using Bowtie2 and BWA
     * Alignment files should be in the alignment folder
     */
    public void mkFastq () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/alignment/small.fq";
        new FeiUtils ().convertTagCount2Fastq(inputFileS, outputFileS);
    }
    
    public void mkSmallTagCount () {
        String inputFileS = "M:/GBStest/tagCount/434GFAAXX_s_4.cnt";
        String outputFileS = "M:/GBStest/tagCount/small.cnt";
        new FeiUtils ().mkSmallTagCountsFile(inputFileS, outputFileS, 10001, 500);
    }
    
    public void pipelinePE () {
        //parseQseq();
        //this.parseFastq();
        //this.checkPETagCounts(); //for checking, not in the pipeline.
        //this.mergePETagCounts();
        this.contigPETagCounts();
        //this.mkPEstatistics();//for presentation, not included in pipeline.
        //this.alignmentStep1();
        //this.alignmentStep2();
        //this.checkAlignmentOfPE(); //check alignment of longer sequence
    }
    
    public void checkAlignmentOfPE () {
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        String topmFileS = "N:/Zea/AllZeaBuild_2.X/04_TOPM/2.6_production/02_MergedTOPM/AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130425.topm";
        String ptopmFileS = "M:/af/ptopm/merge.ptopm";
        String compareTableS = "E:/Research/af/alignmentImprovement/alignmentCompare.txt";
        FeiUtils fu = new FeiUtils ();
        fu.mkAlignmentCompareTable(TOGMFileS, topmFileS, ptopmFileS, compareTableS);
    }
    
    public void alignmentStep2 () {
        String fSamFileS = "M:/af/alignment/f.sam";
        String bSamFileS = "M:/af/alignment/b.sam";
        String contigSamFileS = "M:/af/alignment/c.sam";
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        PETagsOnPhysicalMap topm = new PETagsOnPhysicalMap (ptc, fSamFileS, bSamFileS, contigSamFileS);
        String ptopmFileS = "M:/af/ptopm/merge.ptopm";
        topm.writeDistFile(ptopmFileS, FilePacking.Text);
    }
    
    public void alignmentStep1 () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        String fastaFFileS = "M:/af/alignment/f.fasta.txt";
        String fastaBFileS = "M:/af/alignment/b.fasta.txt";
        String fastaCFileS = "M:/af/alignment/c.fasta.txt";
        ptc.mkFastaFile(fastaFFileS, fastaBFileS, fastaCFileS);
    }
    
    public void mkPEstatistics () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        System.out.println("Tag count is " + ptc.getTagCount());
        System.out.println("Read count is " + ptc.getTotalReadCount());
        System.out.println("Contig count is " + ptc.getContigCount());
        String staFileS = "E:/Research/af/PEstatistics/sta.txt";
        int maxCount = 1;
        for (int i = 0; i < ptc.getTagCount(); i++) {
            if (ptc.getReadCount(i) > maxCount) maxCount = ptc.getReadCount(i);
            //System.out.println(maxCount);            
        }
        System.out.println("Max read count is " + maxCount);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(staFileS), 65536);
            bw.write("FLength\tBLength\tContigLength\tReadCount");
            bw.newLine();
            for (int i = 0; i < 100000; i++) {
                //if (ptc.getContigLengthInLong(i) == 0) continue; 
                bw.write(String.valueOf(ptc.getTagFLength(i))+"\t"+String.valueOf(ptc.getTagBLength(i))+"\t"+String.valueOf(ptc.getContigLength(i))+"\t"+String.valueOf(ptc.getReadCount(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void contigPETagCounts () {
        String infileS = "M:/af/mergePETagCounts/merge.pe.cnt";
        String outfileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        String arguments = "-i " + infileS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        ContigPETagCountPlugin m = new ContigPETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null);
    }
    
    public void mergePETagCounts () {
        String inputDirS = "M:/af/PETagCounts/";
        String outfileS = "M:/af/mergePETagCounts/merge.pe.cnt";
        String arguments = "-i " + inputDirS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        MergePETagCountPlugin m = new MergePETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null); 
    }

    public void parseFastq () {
        String infile1 = "M:/af/Illumina/ImputationP15_1_1_fastq.txt";
        String infile2 = "M:/af/Illumina/ImputationP15_1_2_fastq.txt";
        String keyfile = "M:/af/key/ImputationP15_key.txt";
        String outputDirS = "M:/af/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -l 8 -o " + outputDirS;
        String[] args = arguments.split(" ");
        FastqToPETagCountPlugin q = new FastqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
    public void checkPETagCounts () {
        String outputDirS = "M:/af/PETagCounts/";
        String txtDirS = "M:/af/text/";
        File[] files = new File (outputDirS).listFiles();
        for (int i = 0; i < files.length; i++) {
            PETagCounts p = new PETagCounts (files[i].getAbsolutePath(), FilePacking.Bit);
            String out = txtDirS + "/" + files[i].getName();
            p.writeDistFile(out, FilePacking.Text, 0);
        }
    }
     
    public void parseQseq () {
        String infile1 = "M:/af/Illumina/81546ABXX_8_1_qseq.txt";
        String infile2 = "M:/af/Illumina/81546ABXX_8_2_qseq.txt";
        String keyfile = "M:/af/key/81546ABXX_key.txt";
        String outputDirS = "M:/af/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -o " + outputDirS;
        String[] args = arguments.split(" ");
        QseqToPETagCountPlugin q = new QseqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
}
