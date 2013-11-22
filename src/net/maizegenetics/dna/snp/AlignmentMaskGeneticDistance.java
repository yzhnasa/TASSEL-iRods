/*
 * AlignmentMaskGeneticDistance
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.popgen.distance.IBSDistanceMatrix;
import net.maizegenetics.taxa.Taxon;

import java.awt.*;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author terry
 */
public class AlignmentMaskGeneticDistance extends AbstractAlignmentMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private Map<Integer, Byte> myCache = new LinkedHashMap<Integer, Byte>() {
        protected boolean removeEldestEntry(Map.Entry<Integer, Byte> eldest) {
            return size() > 100;
        }
    };
    private final int myTaxonReference;
    private Alignment myTBitAlignment = null;

    private AlignmentMaskGeneticDistance(Alignment align, int taxonReference, String name, Color color) {
        super(align, name, color, AlignmentMask.MaskType.reference);
        myTaxonReference = taxonReference;
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, Taxon id) {
        List<Integer> index = align.getTaxaList().getIndicesMatchingTaxon(id);
        if ((index == null) || (index.size() == 0)) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index.get(0));
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, String id) {
        List<Integer> index = align.getTaxaList().getIndicesMatchingTaxon(id);
        if ((index == null) || (index.size() == 0)) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index.get(0));
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, int index) {
        if ((index < 0) || (index >= align.getSequenceCount())) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown index: " + index);
        }
        String name = align.getTaxaName(index) + " Genetic Distance";
        return new AlignmentMaskGeneticDistance(align, index, name, null);
    }

    @Override
    public byte getMask(int taxon, int site) {

        Byte result = myCache.get(taxon);
        if (result != null) {
            return result;
        }

//        if (myTBitAlignment == null) {
//            myTBitAlignment = AlignmentUtils.optimizeForTaxa(myAlignment);
//        }

        result = (byte) (IBSDistanceMatrix.computeHetBitDistances(myTBitAlignment, taxon, myTaxonReference)[0] * 255.0);
        myCache.put(taxon, result);
        return result;

    }
}
