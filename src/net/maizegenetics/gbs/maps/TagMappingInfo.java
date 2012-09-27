/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

/**
 *
 * @author edbuckler
 */
public class TagMappingInfo {
    int chromosome=Integer.MIN_VALUE;  // 4 bytes
    byte strand=Byte.MIN_VALUE; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    int startPosition=Integer.MIN_VALUE;  // chromosomal position of the barcoded end of the tag  // 4 bytes
    int endPosition=Integer.MIN_VALUE;  // chromosomal position of the common adapter end of the tag (smaller than startPosition if tag matches minus strand)  // 4 bytes
    byte divergence=Byte.MIN_VALUE;  // number of diverging bp (edit distance) from reference, unknown = Byte.MIN_VALUE
    byte dcoP=Byte.MIN_VALUE;
    byte mapP=Byte.MIN_VALUE;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    
    public TagMappingInfo() {   
    }
    
    public TagMappingInfo(int chromosome, byte strand, int startPosition, 
                int endPosition, byte divergence) {
        this.chromosome=chromosome;
        this.strand=strand;
        this.startPosition=startPosition;
        this.endPosition=endPosition;
        this.divergence=divergence;

    } 
    
    public TagMappingInfo(int chromosome, byte strand, int startPosition, 
                int endPosition, byte divergence, byte[] variantPosOff, byte[] variantDef,
                byte dcoP, byte mapP) {
        this(chromosome, strand, startPosition, endPosition, divergence);
        this.dcoP=dcoP;
        this.mapP=mapP;
    }
        

}
