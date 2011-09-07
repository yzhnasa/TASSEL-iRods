package net.maizegenetics.baseplugins.alignment;


import java.io.Serializable;

import net.maizegenetics.prefs.TasselPrefs;


/**
 * A sequence.
 *
 * @author ajf25 */
public class Sequence implements Serializable {
    // sequence start index
    private int startindex;
    // name
    private String name;
    // sequence
    private String sequence;
    // matched sequence -> can have "."'s
    private String matchedSequence;
    // color sequence based on top sequence
    private String colorTopSequence;
    // whether or not to use the color sequence based on top sequence
    private boolean useColorTopSequence;
    // color sequence based on consensus
    private String colorConsensus;
    // whether or not to use the colored sequence based on consensus
    private boolean useColorConsensus;
    // whether or not to use the matched sequence for displaying purposes
    private boolean useMatchedSequence;
    
    // a quality score for each site in the sequence
    private int[] qualityScore;
    
    /**
     * Constructors
     */
    public Sequence(String name) {
        this(name, "");
    }
    
    public Sequence(String name, String sequence) {
        this.name = name.trim();
        this.sequence = sequence.trim();
        this.matchedSequence = "";
        this.useMatchedSequence = false;
        this.colorTopSequence = "";
        this.useColorTopSequence = false;
        this.colorConsensus = "";
        this.useColorConsensus = false;
        this.startindex = 1;
    }
    
    public Sequence(String name, String sequence, int[] qualityScoreIn) {
        this(name, sequence);
        this.qualityScore = qualityScoreIn;
        
    }
    public int length() {
        return this.sequence.length();
    }
    
    /**
     * Accessors
     */
    // name
    public String getName() {
        return name;
    }
    
    // get the annotation sequence
    public String getSequence() {
        return sequence;
    }
    
    // get the matched annotation sequence
    public String getMatchedSequence() {
        return this.matchedSequence;
    }
    
    public String getColorTopSequence() {
        return this.colorTopSequence;
    }
    
    public String getColorConsensus() {
        return this.colorConsensus;
    }
    
    // get the aligned index of a raw index
    public int getStartIndex() {
        return startindex;
    }
    
    public void setStartIndex(int newStartIndex) {
        this.startindex = newStartIndex;
    }
    
    // set the name
    public void setName(String name) throws Exception {
        this.name = name;
    }
    
    // set the annotation sequence
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    
    // set the matched annotation sequence
    public void setMatchedSequence(String matchedSequence) {
        this.matchedSequence = matchedSequence;
    }
    
    public void setColorTopSequence(String seq) {
        this.colorTopSequence = seq;
    }
    
    public void setColorConsensus(String con) {
        this.colorConsensus = con;
    }
    
    public void setUseQualityScore(boolean use){
        TasselPrefs.putAlignPluginShowQualscore(use);
    }
    
    // set whether or not to display the matched sequence
    public void useMatchedSequence(boolean use) {
        this.useMatchedSequence = use;
    }
    
    public void useColorTopSequence(boolean use) {
        this.useColorTopSequence = use;
    }
    
    public void useColorConsensus(boolean use) {
        this.useColorConsensus = use;
    }
    
    public boolean displayMatchedSequence() {
        return this.useMatchedSequence;
    }
    
    public boolean displayColorTopSequence() {
        return this.useColorTopSequence;
    }
    
    public boolean displayColorConsensus() {
        return this.useColorConsensus;
    }
    
    public boolean displayQualityScore() {
        return TasselPrefs.getAlignPluginShowQualscore();
    }
    
    public int getQualityScore(int i){
        if(i < 0 || i > qualityScore.length - 1 || qualityScore == null){
            return -1;
        }
        return qualityScore[i];
    }
    
    public int[] getQualityScores(){
        return qualityScore;
    }
    
    public void setQualityScores(int[] qualScores){
        qualityScore = qualScores;
    }
}
