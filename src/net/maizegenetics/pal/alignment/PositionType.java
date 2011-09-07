/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.pal.alignment;

/**
 * This class will subset annotated alignment based on a wide range of categories
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class PositionType {

    /** Position Type Group*/
    final static public byte ALL_GROUP = 0;
    final static public byte SILENT_GROUP = 1;
    final static public byte SYNONYMOUS_GROUP = 2;
    final static public byte NONCODING_GROUP = 3;
    final static public byte NONTRANSSCRIBED_GROUP = 4;
    final static public byte INTRON_GROUP = 5;
    final static public byte INDEL_GROUP = 6;
    final static public byte NONCODINGINDEL_GROUP = 7;
    final static public byte NONSYNONYMOUS_GROUP = 8;
    final static public byte CODING_GROUP = 9;
    final static public byte CODINGINDEL_GROUP = 10;
    final static public byte TRANSCRIBED_GROUP = 11;
    final static public String[] GroupText = {"All", "Silent", "Synonymous", "Noncoding",
        "Nontranscribed", "Intron", "Indel", "Noncoding Indel", "Nonsynonymous", "Coding",
        "Coding Indel", "Transcribed"};
    /** Position Type*/
    final static public byte NONTRANSCRIBED_TYPE = 'N';
    /** Position Type*/
    final static public byte ANON_CODING_TYPE = 'C';
    /** Position Type*/
    final static public byte CODON1_TYPE = '1';
    /** Position Type*/
    final static public byte CODON2_TYPE = '2';
    /** Position Type*/
    final static public byte CODON3_TYPE = '3';
    /** Position Type*/
    final static public byte INTRON_TYPE = 'I';

    
    private PositionType() {
    }

    /** Returns a stripped aligment based on the position type
     *  @param positionType the type of position (for example SILENT or SYNONYMOUS)
     *  @return a stripped annotation alignment
     *  */
    public static final Alignment getStrippedAlignment(Alignment aa, int positionType) {
        switch (positionType) {
            case ALL_GROUP: {
                return aa;
            }
            /*      case :
            break;
            default:
            break;*/
        }

        return null;
    }

    /** Returns an array of position types to be included
     *  @param positionType the type of position (for example SILENT or SYNONYMOUS)
     *  */
    public static final char[] getIncludedType(int positionType) {
        switch (positionType) {
            case INDEL_GROUP:
            case ALL_GROUP: {
                char[] inc = {'\u0000', NONTRANSCRIBED_TYPE, CODON1_TYPE, CODON2_TYPE, CODON3_TYPE, INTRON_TYPE, ANON_CODING_TYPE};
                return inc;
            }
            case TRANSCRIBED_GROUP: {
                char[] inc = {CODON1_TYPE, CODON2_TYPE, CODON3_TYPE, INTRON_TYPE, ANON_CODING_TYPE};
                return inc;
            }
            case NONCODINGINDEL_GROUP:
            case NONCODING_GROUP: {
                char[] inc = {NONTRANSCRIBED_TYPE, INTRON_TYPE};
                return inc;
            }
            case INTRON_GROUP: {
                char[] inc = {INTRON_TYPE};
                return inc;
            }
            case SILENT_GROUP: {
                char[] inc = {NONTRANSCRIBED_TYPE, CODON3_TYPE, INTRON_TYPE};
                return inc;
            }
            case SYNONYMOUS_GROUP:
            case NONSYNONYMOUS_GROUP:
            case CODINGINDEL_GROUP:
            case CODING_GROUP: {
                char[] inc = {CODON1_TYPE, CODON2_TYPE, CODON3_TYPE, ANON_CODING_TYPE};
                return inc;
            }
            case NONTRANSSCRIBED_GROUP: {
                char[] inc = {NONTRANSCRIBED_TYPE};
                return inc;
            }
            default: {
                return null;
            }
        }
    }

    public static final boolean isExon(byte t) {
        if ((t == CODON1_TYPE) || (t == CODON2_TYPE) || (t == CODON3_TYPE) || (t == ANON_CODING_TYPE)) {
            return true;
        }
        return false;
    }

    /** Returns an array of position types to be included
     *  @param positionType the type of position (for example SILENT or SYNONYMOUS)
     *  */
    public static final String getPositionGroupName(int positionType) {
        if (positionType < GroupText.length) {
            return GroupText[positionType];
        } else {
            return "Group Not Found";
        }
    }
}