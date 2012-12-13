/*
 * AlignmentMaskGeneticDistance
 */
package net.maizegenetics.pal.alignment;

import java.awt.Color;
import java.util.WeakHashMap;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry
 */
public class AlignmentMaskGeneticDistance extends AbstractAlignmentMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private WeakHashMap<Integer, Byte> myCachedDistances = new WeakHashMap<Integer, Byte>(100);
    private final int myTaxonReference;
    private final Alignment myAlignment;

    private AlignmentMaskGeneticDistance(Alignment align, int taxonReference, String name, Color color) {
        super(align, name, color, AlignmentMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, Identifier id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, String id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
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

        Byte result = myCachedDistances.get(taxon);
        if (result != null) {
            return result;
        }

        result = (byte) (IBSDistanceMatrix.computeHetBitDistances(myAlignment, taxon, myTaxonReference) * 255.0);
        myCachedDistances.put(taxon, result);
        return result;

    }
}
