package net.maizegenetics.baseplugins.alignment;

//import jalview.AlignmentPanel;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.SimpleAlignment;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.datatype.TextDataType;

import javax.swing.*;
import java.awt.*;


/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jun 3, 2004
 * Time: 10:27:31 AM
 */
public class AlignmentFrame extends JPanel {
    // sequence stuff
    // store the annotation alignment so changes can be pushed back

    Alignment aa;
    Sequence[] sequences;
    StringBuffer alignedWithHeadersSB;
    int maxSeqLength;
    ScrollableViewport alignment_name_vp, alignment_col_vp, alignment_data_vp;
    AlignmentNamePanel alignment_name_panel;
    AlignmentGroupPanel alignment_col_panel;
    AlignmentPanel alignment_data_panel;
    public JPanel alignment_overall_panel;
    JScrollBar alignment_vert_bar, alignment_hor_bar;
    JFrame theParent;

    // offset to handle the different number of blank JLabels() placed in the namePanel
    // in order to align the sequence name and sequence correctly
    int nameComponentOffset;
    FontMetrics theFM;
    public int headerOffset = 40;
    public int xOffset = 7;   //the width of a character
    public int yOffset = 16;  //the height of a character

    public AlignmentFrame(Alignment aa) {
        super();
        this.aa = aa;
        alignedWithHeadersSB = null;
        nameComponentOffset = 2;
        int seqCount = aa.getIdGroup().getIdCount();
        sequences = new Sequence[seqCount];
        for (int i = 0; i < seqCount; i++) {
            sequences[i] = new Sequence(aa.getIdGroup().getIdentifier(i).toString(), aa.getAlignedSequenceString(i).trim());
        }
        init(null);
    }

    // this is one of the mainly used constructors
    public AlignmentFrame(Alignment aa, JFrame theParent) {
        super();
        this.aa = aa;
        this.theParent = theParent;
        theFM = theParent.getGraphics().getFontMetrics(new Font("Monospaced", Font.PLAIN, 12));
        xOffset = theFM.charWidth('G');
        yOffset = theFM.getHeight();
        headerOffset = yOffset * 3;
        System.out.println("theFM.getWidths()=" + xOffset + " theFM.getHeight()=" + yOffset);
        alignedWithHeadersSB = null;
        nameComponentOffset = 2;
        int seqCount = aa.getIdGroup().getIdCount();
        sequences = new Sequence[seqCount];
        for (int i = 0; i < seqCount; i++) {
            int siteCount = aa.getSiteCount();
            int[] qualScore = new int[siteCount];
            //if (aa.hasQualityScores()) {
            //    for (int site = 0; site < siteCount; site++) {
            //        qualScore[site] = aa.getQualityScore(i, site);
            //    }
            //}
            //      sequences[i] = new Sequence( aa.getIdentifier(i).toString(), aa.getAlignedSequenceString(i).trim() , qualScore);
            sequences[i] = new Sequence(aa.getIdGroup().getIdentifier(i).toString(), aa.getFormattedSequenceString(i).trim(), qualScore);
        }
        init(null);
    }

    // this constructor is used for when alignWithHeaders would be called over
    // quickAlignWithHeaders in the AlignmentPrintUtils method alignWithHeadersToStringBuffer
    public AlignmentFrame(Alignment aa, StringBuffer sb, JFrame theParent) {
        super();
        this.aa = aa;
        this.theParent = theParent;
        theFM = theParent.getGraphics().getFontMetrics(new Font("Monospaced", Font.PLAIN, 12));
        xOffset = theFM.charWidth('G');
        yOffset = theFM.getHeight();
        headerOffset = yOffset * 5;
        alignedWithHeadersSB = sb;
        nameComponentOffset = 5;
        int seqCount = aa.getIdGroup().getIdCount();
        sequences = new Sequence[seqCount];
        // parse the sb for the sequences
        // format of sb is:
        // pos>
        // pos>
        // pos>
        // pos>
        // pos>
        // Taxa blah blah \n
        // name    seq
        int index = sb.indexOf("Taxa");
        index = sb.indexOf("\n", index);
        System.out.println("index = " + index);



        // while we have not processed all the names and seqs
//        for(int i = 0; i < seqCount; i++)
//        {
//            int indexSpace = sb.indexOf(" ", index+1);
//            String testString = sb.substring(0, index+1).toString();
//            System.out.println("testString = " + testString);
//            System.out.println("indexSpace of " + i + " = " + indexSpace);
//            int indexLine = -1;
//            String name = null;
//            String seq = null;
//            if(indexSpace == -1){      // this protects against  throwing a StringIndexOutOfBoundsException
//                                       //  as  a result of variations in formatting
//                name = "error";
//                seq  = "GCTA";
//            }else{
//                name = sb.substring(index+1, indexSpace);
//                System.out.println("name = " + name);
//                //todo at times is unable to parse this & throws java.lang.StringIndexOutOfBoundsException
//                indexLine = sb.indexOf("\n", indexSpace);
//                seq = sb.substring(indexSpace, indexLine).trim();
//            }
//
//            sequences[i] = new Sequence(name, seq);
//            index = indexLine;
//        }

        // This is a temporary replacement of the above code.  The above problem actually has to do with
        // the fact that AJ never was able to get the changes in the AlignmentFrame to be stored on the
        // data tree.  He unsuccessfully attempted to parse the results.  This unsuccessful code path was
        // never used until the Taxa filter was started to work with sequence data.  For the time being,
        // it is more important to be able to subset taxa than to be able to do minor tweaking on alignments. -DEK
        //get the names and sequences from the aa
        int seqLength = aa.getSiteCount();
        for (int i = 0; i < seqCount; i++) {
            String name = aa.getIdGroup().getIdentifier(i).getName();
            char[] seqArray = new char[seqLength];
            for (int j = 0; j < seqLength; j++) {
                seqArray[j] = aa.getBaseChar(i, j);
            }
            String sequence = new String(seqArray);
            sequences[i] = new Sequence(name, sequence);
        }
        init(sb);
    }

    public AlignmentFrame() {
        super();
        this.aa = null;
        //**************************************************************************************
        // for testing purposes - create the sequences
        sequences = new Sequence[100];
        String test = "ACTGACTGACTGACTGACTGAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTACTGACTGACTGACTGACTGAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTACTGACTGACTGACTGACTGAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTACTGACTGACTGACTGACTGAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTEND";
        for (int i = 0; i < 100; i++) {
            sequences[i] = new Sequence("hello world " + (i + 1), test);
        }  // for testing purposes
        //**************************************************************************************
        init(null);
    }

    public void init(StringBuffer sb) {
        recalculateMaxSeqLength();

        // init everything
        alignment_name_vp = new ScrollableViewport();
        alignment_col_vp = new ScrollableViewport();
        alignment_data_vp = new ScrollableViewport();
        if (sb == null) {
            alignment_name_panel = new AlignmentNamePanel(this, sequences, theFM);
            alignment_col_panel = new AlignmentGroupPanel(this, sequences, theFM);
        } else {
            alignment_name_panel = new AlignmentNamePanel(this, sequences, sb, theFM);
            alignment_col_panel = new AlignmentGroupPanel(this, sequences, sb, aa, theFM);
        }
        alignment_data_panel = new AlignmentPanel(this, sequences, theFM);
        alignment_overall_panel = new JPanel();
        alignment_vert_bar = new JScrollBar(JScrollBar.VERTICAL);
        alignment_hor_bar = new JScrollBar(JScrollBar.HORIZONTAL);
        //System.out.println("maxseqlength before scrollbar created is " + maxSeqLength);    //****
        //alignment_hor_bar = new JScrollBar(JScrollBar.HORIZONTAL, 0, 85, 0, maxSeqLength );

        // set up name viewport
        alignment_name_vp.putClientProperty("EnableWindowBlit", null);
        alignment_name_vp.setBackground(Color.white);
        alignment_name_vp.setVerticalScrollbar(alignment_vert_bar);
        alignment_name_vp.setView(alignment_name_panel);

        // set up col viewport
        alignment_col_vp.putClientProperty("EnableWindowBlit", null);
        alignment_col_vp.setBackground(Color.white);
        alignment_col_vp.setHorizontalScrollbar(alignment_hor_bar);
        alignment_col_vp.setView(alignment_col_panel);

        // set up data viewport
        alignment_data_vp.putClientProperty("EnableWindowBlit", null);
        alignment_data_vp.setBackground(Color.white);
        alignment_data_vp.setVerticalScrollbar(alignment_vert_bar);
        alignment_data_vp.setHorizontalScrollbar(alignment_hor_bar);
        alignment_data_vp.setView(alignment_data_panel);

        // create a panel to handle the layout in the CENTER region of the overall_panel
        JPanel center_panel = new JPanel();
        center_panel.setLayout(new BorderLayout());

        // add all the CENTER components to the correct region in center_panel
        center_panel.add(BorderLayout.NORTH, alignment_col_vp);
        center_panel.add(BorderLayout.CENTER, alignment_data_vp);

        // set the layout for the alignment_overall_panel and add all the viewports and scrollbars
        alignment_overall_panel.setLayout(new BorderLayout());
        alignment_overall_panel.add(BorderLayout.WEST, alignment_name_vp);
        alignment_overall_panel.add(BorderLayout.CENTER, center_panel);
        alignment_overall_panel.add(BorderLayout.EAST, alignment_vert_bar);
        alignment_overall_panel.add(BorderLayout.SOUTH, alignment_hor_bar);

        add(alignment_overall_panel);

        updateAnnotationAlignment();
    }

    public void updateAnnotationAlignment() {
        // don't allow saving when viewing haplotypes and amino acids
        if (alignedWithHeadersSB != null) {
            return;
        }

        // create a String[] of the sequence names, and sequences themselves
        String[] names = new String[sequences.length];
        String[] seq = new String[sequences.length];
        float[][] siteScores = new float[sequences.length][];

        for (int i = 0; i < sequences.length; i++) {
            names[i] = sequences[i].getName();
            seq[i] = sequences[i].getSequence();
            siteScores[i] = new float[seq[i].length()];
            for (int j = 0; j < seq[i].length(); j++) {
                siteScores[i][j] = sequences[i].getQualityScore(j);
            }
        }
        SimpleIdGroup theSimpleIdGroup = new SimpleIdGroup(names);
        // this is creating a new annotationalignment and updating the datum stored in the book
        // located in the DataTreePanel
        //this.aa = new SimpleAlignment(theSimpleIdGroup, seq, qualScores, new Nucleotides());
        this.aa = new SimpleAlignment(theSimpleIdGroup, seq, new TextDataType(), null, null, null, aa.getLocusName(0), siteScores, null);
    //Ed is turning off this updateDatum as we are loosing locus name information, and perhaps all annotations
    //this.book.updateDatum(this.aa);
    }

    public void swapSequences(int start, int end) {
        //System.out.println("Start and End are " + start + " and " + end);    //****
        //System.out.println("Swapping indices of " + (start-1) + " and " + (end-1));    //****

        // grab the 2 sequencecomponents
        SequenceComponent scStart = (SequenceComponent) alignment_data_panel.getComponent(start);
        SequenceComponent scEnd = (SequenceComponent) alignment_data_panel.getComponent(end);

        // swap the 2 sequencecomponents sequences
        Sequence scStartSeq = scStart.getSequence();
        Sequence scEndSeq = scEnd.getSequence();
        scStart.setSequence(scEndSeq);
        scEnd.setSequence(scStartSeq);

        // swap the actual sequences
        Sequence temp = sequences[start - 1];
        sequences[start - 1] = sequences[end - 1];
        sequences[end - 1] = temp;

        // swap the 2 name labels
        SequenceNameComponent sncStart = (SequenceNameComponent) alignment_name_panel.getComponent(start + nameComponentOffset);
        SequenceNameComponent sncEnd = (SequenceNameComponent) alignment_name_panel.getComponent(end + nameComponentOffset);
        Sequence sncStartSeq = sncStart.getSequence();
        Sequence sncEndSeq = sncEnd.getSequence();
        //System.out.println("Starting name is " + sncStartSeq.getName() );    //****
        //System.out.println("Ending name is " + sncEndSeq.getName() );    //****
        sncStart.setSequence(sncEndSeq);
        sncEnd.setSequence(sncStartSeq);
    }
    /*
    private String removeEndingDashes(String seq)
    {
    StringBuffer sb = new StringBuffer(seq);

    while( sb.length() > 0 && sb.charAt(sb.length()-1) == Alignment.GAP )
    sb.deleteCharAt(sb.length()-1);
    return sb.toString();
    }
     */

    private void recalculateMaxSeqLength() {
        int maxLen = -1;
        for (int i = 0; i < sequences.length; i++) {
            Sequence tempSeq = sequences[i];
            int tempLen = tempSeq.length();
            if (tempLen > maxLen) {
                maxLen = tempLen;
            }
        }
        maxSeqLength = Math.max(maxLen, 10);
    }
}
