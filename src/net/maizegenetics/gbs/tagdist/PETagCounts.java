/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.tagdist;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class PETagCounts extends AbstractPETags {
    int[] readCount;
    
    public PETagCounts (String inFile, FilePacking format) {
        this.readDistFile(inFile, format);
    }
    
    public PETagCounts (int tagLengthInLong, int tagNum,String[] taxaName) {
        this.iniMatrix(tagLengthInLong, tagNum, taxaName.length);
        this.taxaName = taxaName;
    }
    
    public PETagCounts getCollapsedPETagCounts () {
        int tagNum = this.getTagCount() - this.collapseCounts();
        PETagCounts petc = new PETagCounts (this.getTagSizeInLong(), tagNum, this.taxaName);
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.getReadCount(i) == 0) continue;
            for (int j = 0; j < this.getTagSizeInLong(); j++) {
                petc.tagsF[j][cnt] = this.tagsF[j][i];
                petc.tagsB[j][cnt] = this.tagsB[j][i];
            }
            petc.tagFLength[cnt] = this.tagFLength[i];
            petc.tagBLength[cnt] = this.tagBLength[i];
            petc.readCount[cnt] = this.readCount[i];
            cnt++;
        }
        return petc;
    }
    
    public int collapseCounts () {
        this.sort();
        int collapsedRows = 0;
        for (int i = 1; i < this.getTagCount(); i++) {
            if (this.compare(i, i-1) == 0) {
                readCount[i] += readCount[i - 1];
                readCount[i - 1] = 0;
                collapsedRows++;
            }
        }
        System.out.println("Rows collapsed:" + collapsedRows);
        return collapsedRows;
    }
    
    public PETagCounts getMergedPETagCounts (PETagCounts another, boolean ifCollapsed) {
        if (!ifCollapsed) another = another.getCollapsedPETagCounts();
        TreeSet<String> nameSet = new TreeSet();
        for (int i = 0; i < this.getTaxaCount(); i++) nameSet.add(this.getTaxaName(i));
        for (int i = 0; i < another.getTaxaCount(); i++) nameSet.add(another.getTaxaName(i));
        String[] taxaName = nameSet.toArray(new String[nameSet.size()]);
        PETagCounts petc = new PETagCounts(this.tagLengthInLong, this.getTagCount()+another.getTagCount(), taxaName);
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            for (int j = 0; j < this.getTagSizeInLong(); j++) {
                petc.tagsF[j][cnt] = this.tagsF[j][i];
                petc.tagsB[j][cnt] = this.tagsB[j][i];
            }
            petc.tagFLength[cnt] = this.tagFLength[i];
            petc.tagBLength[cnt] = this.tagBLength[i];
            petc.readCount[cnt] = this.readCount[i];
            cnt++;
        }
        for (int i = 0; i < another.getTagCount(); i++) {
            for (int j = 0; j < another.getTagSizeInLong(); j++) {
                petc.tagsF[j][cnt] = another.tagsF[j][i];
                petc.tagsB[j][cnt] = another.tagsB[j][i];
            }
            petc.tagFLength[cnt] = another.tagFLength[i];
            petc.tagBLength[cnt] = another.tagBLength[i];
            petc.readCount[cnt] = another.readCount[i];
            cnt++;
        }
        return petc.getCollapsedPETagCounts();
    }
    
    public int getReadCount (int index) {
        return readCount[index];
    }
    
    private int getTagNumWithMincount (int minCount) {
        int num = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (readCount[i] >= minCount) num++;
        }
        return num;
    }
    
    @Override
    public void readDistFile (String infileS, FilePacking format) {
        System.out.println("Reading PETagCounts file to " + infileS);
        File infile = new File (infileS);
        switch (format) {
            case Text:
                readTextPETagCountsFile(infile);
                break;
            default:
                readBinaryPETagCountsFile(infile);
                break;
        }
        System.out.println("PETagCounts file read");
    }
    
    private void readBinaryPETagCountsFile (File infile) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            int taxaNum = dis.readInt();
            tagLengthInLong = dis.readInt();
            int currentSize = 8;
            for (int i = 0; i < taxaNum; i++) {
                currentSize = dis.readUTF().length();
            }
            int lineSize = (tagLengthInLong*8+2)*2+4;
            int tagNum = (int)((infile.length()-currentSize)/lineSize);
            dis.close();
            dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            taxaNum = dis.readInt();
            tagLengthInLong = dis.readInt();
            this.iniMatrix(tagLengthInLong, tagNum, taxaNum);
            for (int i = 0; i < taxaNum; i++) {
                taxaName[i] = dis.readUTF();
            }
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsF[j][i] = dis.readLong();
                }
                tagFLength[i] = dis.readShort();
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsB[j][i] = dis.readLong();
                }
                tagBLength[i] = dis.readShort();
                readCount[i] = dis.readInt();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void readTextPETagCountsFile (File infile) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(infile), 65536);
            int taxaNum = Integer.valueOf(br.readLine());
            tagLengthInLong = Integer.valueOf(br.readLine());
            int tagNum = Integer.valueOf(br.readLine());
            this.iniMatrix(tagLengthInLong, tagNum, taxaNum);
            taxaName = br.readLine().split("\t");
            for (int i = 0; i < tagNum; i++) {
                String[] temp = br.readLine().split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(temp[0]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsF[j][i] = t[j];
                }
                tagFLength[i] = Byte.valueOf(temp[1]);
                t = BaseEncoder.getLongArrayFromSeq(temp[2]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsB[j][i] = t[j];
                }
                tagFLength[i] = Byte.valueOf(temp[3]);
                readCount[i] = Integer.valueOf(temp[4]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void writeDistFile (String outfileS, FilePacking format, int minCount) {
        System.out.println("Writing PETagCounts file to " + outfileS);
        int outTagNum = this.getTagNumWithMincount(minCount);
        switch (format) {
            case Text:
                writeTextPETagCountsFile(outfileS, outTagNum, minCount);
                break;
            default:
                writeBinaryPETagCountsFile(outfileS, minCount);
                break;
        }
        System.out.println("PETagCounts file written");
    }
    
    private void writeBinaryPETagCountsFile (String outfileS, int minCount) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(this.getTaxaCount());
            dos.writeInt(tagLengthInLong);
            for (int i = 0; i < this.getTaxaCount(); i++) {
                dos.writeUTF(this.taxaName[i]);
            }
            for (int i = 0; i < this.getTagCount(); i++) {
                if (readCount[i] < minCount) continue;
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tagsF[j][i]);
                }
                dos.writeShort(this.tagFLength[i]);
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tagsB[j][i]);
                }
                dos.writeShort(this.tagBLength[i]);
                dos.writeInt(this.readCount[i]);
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void writeTextPETagCountsFile (String outfileS, int outTagNum, int minCount) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write(String.valueOf(this.getTaxaCount()));
            bw.newLine();
            bw.write(String.valueOf(tagLengthInLong));
            bw.newLine();
            bw.write(String.valueOf(outTagNum));
            bw.newLine();
            for (int i = 0; i < this.getTaxaCount(); i++) {
                bw.write(this.taxaName[i]+"\t");
            }
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.getReadCount(i) < minCount) continue;
                bw.write(BaseEncoder.getSequenceFromLong(this.getTagF(i))+"\t"+String.valueOf(this.getTagFLength(i))+"\t");
                bw.write(BaseEncoder.getSequenceFromLong(this.getTagB(i))+"\t"+String.valueOf(this.getTagBLength(i))+"\t");
                bw.write(String.valueOf(this.getReadCount(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    @Override
    protected void iniMatrix (int tagLengthInLong, int tagNum, int taxaNum) {
        super.iniMatrix(tagLengthInLong, tagNum, taxaNum);
        readCount = new int[tagNum];
    }
    
    @Override
    public void swap(int index1, int index2) {
         super.swap(index1, index2);
         int temp = readCount[index1];
         readCount[index1] = readCount[index2];
         readCount[index2] = temp;
    }
}
