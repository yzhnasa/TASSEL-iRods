/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;

/**
 *
 * @author Fei Lu
 */
public class FeiPipelines {
    
    public FeiPipelines () {
        //parseQseq();
        this.parseFastq();
        this.checkPETagCounts();
        this.mergePETagCounts();
        this.contigPETagCounts();
    }
    
    public static void main (String[] args) {
        new FeiPipelines();
        
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
        String infile1 = "M:/af/Illumina/test/ImputationP15_1_1_fastq.txt";
        String infile2 = "M:/af/Illumina/test/ImputationP15_1_2_fastq.txt";
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
