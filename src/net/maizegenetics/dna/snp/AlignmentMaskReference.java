/*
 * AlignmentMaskReference
 */
package net.maizegenetics.dna.snp;

import java.awt.Color;
import java.util.List;
import net.maizegenetics.taxa.Taxon;

/**
 *
 * @author terry
 */
public class AlignmentMaskReference extends AbstractAlignmentMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private final int myTaxonReference;
    private final Alignment myAlignment;

    private AlignmentMaskReference(Alignment align, int taxonReference, String name, Color color) {
        super(align, name, color, AlignmentMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align) {
        return getInstanceCompareReference(align, -1);
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align, Taxon id) {
        List<Integer> index = align.taxa().getIndicesMatchingTaxon(id);
        if ((index == null) || (index.size() == 0)) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index.get(0));
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align, String id) {
        List<Integer> index = align.taxa().getIndicesMatchingTaxon(id);
        if ((index == null) || (index.size() == 0)) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index.get(0));
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align, int index) {
        if ((index < -1) || (index >= align.numberOfTaxa())) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown index: " + index);
        }
        String name;
        if (index == -1) {
            name = "Alignment Reference";
        } else {
            name = align.taxaName(index) + " Reference";
        }
        return new AlignmentMaskReference(align, index, name, getNextColor());
    }

    @Override
    public byte getMask(int taxon, int site) {
        if ((myTaxonReference == -1) && (AlignmentUtils.isEqualOrUnknown(myAlignment.genotype(taxon, site), myAlignment.getReferenceAllele(site)))) {
            return 0;
        } else if (AlignmentUtils.isEqualOrUnknown(myAlignment.genotypeArray(taxon, site), myAlignment.genotypeArray(myTaxonReference, site))) {
            return 0;
        } else {
            return 1;
        }
    }
}
