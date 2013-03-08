/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

/**
 *
 * @author edbuckler
 */
public class SiteMappingInfo {
    public int chromosome=Integer.MIN_VALUE;  // 4 bytes
    public byte strand=Byte.MIN_VALUE; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    public int position=Integer.MIN_VALUE;  // chromosomal position of the barcoded end of the tag  // 4 bytes
    public float r2=Float.NaN;
    public float mapP=Float.NaN;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    
    public SiteMappingInfo() {   
    }
    
    public SiteMappingInfo(int chromosome, byte strand, int position, 
                float r2, float mapP) {
        this.chromosome=chromosome;
        this.strand=strand;
        this.position=position;
        this.r2=r2;
        this.mapP=mapP;

    } 
    
        

}
